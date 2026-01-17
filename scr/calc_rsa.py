import os
import sys
import subprocess
import warnings
import re
import time
import tempfile
import logging
import requests
import shutil

# --- 引入依赖库 ---
try:
    from Bio.PDB import PDBParser
    from Bio import BiopythonWarning
    import pandas as pd
    import freesasa
except ImportError as e:
    print(f"环境错误: 缺少必要库 ({e})")
    print("请运行: pip install pandas openpyxl biopython requests freesasa")
    sys.exit(1)

# ================= 用户自定义参数 (CONFIG) =================

# 1. 输入与输出设置
METADATA_FILE = "/Users/zhanghaonan/Documents/ML/plaxy/data/uniprotkb_taxonomy_id_2188_Methanococcus voltae_40c.tsv"
OUTPUT_CSV = "rsa_Methanococcus voltae_40c.tsv"

# --- [新增] 失败 ID 清单文件 (纯 ID，一行一个，方便后续放入预测软件) ---
FAILED_ID_FILE = "failed_tasks_list_Methanococcus voltae_40c.txt"

# --- [保留] 详细错误日志 (包含时间戳和具体原因) ---
ERROR_LOG_FILE = "processing_details_Methanococcus voltae_40c.log"

# 2. Excel/CSV 列名映射
COL_ENTRY = "Entry"          # [必须] Uniprot ID 列
COL_MIN_PH = None        # [可选] 最小 pH 列
COL_MAX_PH = None        # [可选] 最大 pH 列
COL_ORGANISM = None    # [可选] 生物体列
SHEET_NAME = 0               # [仅Excel] Sheet名称或索引

# 3. 筛选条件 (设为 None 则不进行该项筛选)
FILTER_MIN_PH = None
FILTER_MAX_PH = None

# 4. 计算参数
TARGET_PH = 7.0
SLEEP_ON_FAIL = 0.2

# ==========================================================

# 设置日志 (详细版)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(ERROR_LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)

# 屏蔽警告
warnings.simplefilter('ignore', BiopythonWarning)
freesasa.setVerbosity(freesasa.nowarnings)

# 常量定义
MAX_ASA = {
    'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0, 'E': 214.0,
    'Q': 214.0, 'G': 97.0,  'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
    'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0, 'T': 163.0, 'W': 264.0,
    'Y': 255.0, 'V': 165.0
}
PKA_VALUES = {'D': 3.9, 'E': 4.2, 'H': 6.0, 'C': 8.3, 'Y': 10.07, 'K': 10.53, 'R': 12.48}
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'ASH': 'D', 'GLH': 'E', 'LYN': 'K', 'ARN': 'R', 
    'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C', 'CYM': 'C', 'MSE': 'M'
}

def clean_ph_value(val):
    if pd.isna(val): return None
    match = re.search(r"(\d+(\.\d+)?)", str(val).strip())
    return float(match.group(1)) if match else None

def load_metadata(file_path):
    if not os.path.exists(file_path):
        logging.error(f"文件不存在: {file_path}")
        sys.exit(1)
    ext = os.path.splitext(file_path)[1].lower()
    try:
        if ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, sheet_name=SHEET_NAME)
        elif ext == '.csv':
            df = pd.read_csv(file_path)
        elif ext in ['.tsv', '.txt']:
            df = pd.read_csv(file_path, sep='\t')
        else:
            logging.error(f"不支持的文件格式: {ext}")
            sys.exit(1)
        df.columns = df.columns.str.strip()
        return df
    except Exception as e:
        logging.error(f"读取文件失败: {e}")
        sys.exit(1)

def get_primary_accession(old_id):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{old_id}"
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            primary = data.get('primaryAccession', '')
            if primary and primary != old_id:
                return primary
    except Exception:
        pass
    return None

def download_alphafold_structure(uniprot_id):
    uid = uniprot_id.strip()
    # 构造尝试列表：先试原始 ID，再试映射后的主 ID
    versions = ["v6", "v4", "v3"] # 优先试最新的
    
    # 1. 第一阶段：直接尝试下载（最快）
    for v in versions:
        # 尝试常见的两种命名格式
        urls = [
            f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{v}.pdb",
            f"https://alphafold.ebi.ac.uk/files/AF-{uid}-model_{v}.pdb" # 简单编号格式
        ]
        for url in urls:
            try:
                resp = requests.get(url, timeout=5) # 简单下载给 5s 即可
                if resp.status_code == 200:
                    fd, path = tempfile.mkstemp(suffix=".pdb")
                    with os.fdopen(fd, 'wb') as tmp:
                        tmp.write(resp.content)
                    logging.info(f"成功下载: {uid} (版本 {v})")
                    return path, uid
            except: continue

    # 2. 第二阶段：如果失败，再求助 UniProt 映射和官方 API（较慢）
    logging.info(f"直接下载失败，尝试 API 映射: {uid}")
    primary = get_primary_accession(uid)
    search_id = primary if primary else uid
    
    try:
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{search_id}"
        resp = requests.get(api_url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            pdb_url = data[0].get('pdbUrl') if isinstance(data, list) and data else None
            if pdb_url:
                r = requests.get(pdb_url, timeout=15)
                if r.status_code == 200:
                    fd, path = tempfile.mkstemp(suffix=".pdb")
                    with os.fdopen(fd, 'wb') as tmp: tmp.write(r.content)
                    return path, search_id
    except Exception as e:
        logging.debug(f"API 请求异常: {e}")
    
    return None, None

def run_pdb2pqr_temp(pdb_path, ph):
    fd, pqr_path = tempfile.mkstemp(suffix=".pqr")
    os.close(fd)
    cmd = [
        "pdb2pqr", "--ff=AMBER", "--titration-state-method=propka",
        f"--with-ph={ph}", "--noopt", pdb_path, pqr_path
    ]
    try:
        subprocess.run(cmd, capture_output=True, check=True)
        return pqr_path
    except subprocess.CalledProcessError as e:
        logging.warning(f"PDB2PQR 失败: {e.stderr.decode('utf-8')[:200]}...")
        if os.path.exists(pqr_path): os.remove(pqr_path)
        return None

def calculate_hh_charge(res_name_1, ph):
    pka = PKA_VALUES.get(res_name_1)
    if pka is None: return 0.0
    if res_name_1 in ['D', 'E', 'C', 'Y']:
        return -1.0 / (1.0 + 10**(pka - ph))
    elif res_name_1 in ['H', 'K', 'R']:
        return 1.0 / (1.0 + 10**(ph - pka))
    return 0.0

def get_charges_from_pqr(pqr_path):
    charges = {}
    if not pqr_path or not os.path.exists(pqr_path): return {}
    with open(pqr_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                try:
                    charge = float(parts[-2])
                    res_id = None
                    for token in reversed(parts[:-2]):
                        if token.lstrip('-').isdigit():
                            res_id = int(token)
                            break
                    if res_id is not None:
                        charges[res_id] = charges.get(res_id, 0.0) + charge
                except: continue
    return charges

def analyze_structure(pdb_path, pqr_path, target_ph):
    parser = PDBParser(QUIET=True)
    results = []
    try:
        structure = parser.get_structure('tmp', pdb_path)
        f_structure = freesasa.Structure(pdb_path)
        f_result = freesasa.calc(f_structure)
        f_areas = f_result.residueAreas()
        pqr_charges = get_charges_from_pqr(pqr_path)
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                query_chain = chain_id if chain_id != ' ' else 'A'
                for residue in chain:
                    if residue.id[0] != ' ': continue
                    res_id = residue.id[1]
                    res_name_3 = residue.get_resname()
                    res_name_1 = AA_3TO1.get(res_name_3, '?')
                    if res_name_1 == '?': continue
                    sasa = 0.0
                    try:
                        if query_chain in f_areas:
                            sasa = f_areas[query_chain][str(res_id)].total
                        elif 'A' in f_areas:
                            sasa = f_areas['A'][str(res_id)].total
                    except KeyError: pass
                    rsa = min(1.0, sasa / MAX_ASA.get(res_name_1, 1.0))
                    c_pqr = pqr_charges.get(res_id, 0.0)
                    c_avg = calculate_hh_charge(res_name_1, target_ph)
                    results.append({
                        'res_id': res_id, 'res_name': res_name_1,
                        'sasa': sasa, 'rsa': rsa,
                        'charge_pqr': c_pqr, 'charge_avg': c_avg
                    })
    except Exception as e:
        logging.error(f"分析结构时出错: {e}")
        return []
    return results

# --- 新增函数：记录失败ID ---
def record_failed_id(uid):
    """将失败的 ID 追加写入到纯文本清单中"""
    try:
        with open(FAILED_ID_FILE, "a") as f:
            f.write(f"{uid}\n")
    except Exception as e:
        print(f"无法写入失败清单: {e}")

# ================= 主程序 =================

if __name__ == "__main__":
    print("=== 蛋白结构获取与分析自动化脚本 (带失败记录) ===")
    print(f"失败的ID将保存在: {FAILED_ID_FILE}")
    
    df = load_metadata(METADATA_FILE)
    if COL_ENTRY not in df.columns:
        print(f"错误: 找不到必要的列 '{COL_ENTRY}'")
        sys.exit(1)
        
    # ================= 筛选逻辑 (修改版) =================
    # 预处理：生成临时的数值型 pH 列
    if COL_MIN_PH and COL_MIN_PH in df.columns:
        df['temp_min_ph'] = df[COL_MIN_PH].apply(clean_ph_value)
    else:
        df['temp_min_ph'] = pd.NA

    if COL_MAX_PH and COL_MAX_PH in df.columns:
        df['temp_max_ph'] = df[COL_MAX_PH].apply(clean_ph_value)
    else:
        df['temp_max_ph'] = pd.NA

    # 初始化筛选掩码 (Mask)，默认为 False (不保留)
    # 如果用户两个参数都设为 None，则默认保留所有数据
    if FILTER_MIN_PH is None and FILTER_MAX_PH is None:
        mask = pd.Series([True] * len(df), index=df.index)
        print("提示: 未设置任何 pH 筛选条件，保留所有数据。")
    else:
        mask = pd.Series([False] * len(df), index=df.index)
        filter_applied = False

        # 条件 A: min_ph 不得高于 (<=) FILTER_MIN_PH
        if FILTER_MIN_PH is not None and 'temp_min_ph' in df.columns:
            # 注意：排除 NaN 值，只有明确有数值且 <= 阈值的才为 True
            condition_min = (df['temp_min_ph'].notna()) & (df['temp_min_ph'] <= FILTER_MIN_PH)
            mask = mask | condition_min  # 逻辑 OR
            print(f"筛选条件 A 应用: Min pH <= {FILTER_MIN_PH}")
            filter_applied = True

        # 条件 B: max_ph 不得低于 (>=) FILTER_MAX_PH
        # 也就是 "或者" max_ph >= FILTER_MAX_PH
        if FILTER_MAX_PH is not None and 'temp_max_ph' in df.columns:
            condition_max = (df['temp_max_ph'].notna()) & (df['temp_max_ph'] >= FILTER_MAX_PH)
            mask = mask | condition_max  # 逻辑 OR
            print(f"筛选条件 B 应用: Max pH >= {FILTER_MAX_PH}")
            filter_applied = True
        
        # 如果设置了筛选参数但数据列不存在，可能会导致 mask 全 False，这里不做额外处理，直接筛选
        
        # 应用筛选
        original_count = len(df)
        df = df[mask]
        print(f"筛选结果: {original_count} -> {len(df)} (保留满足任一条件的条目)")
    # ====================================================

    filtered_count = len(df)
    print(f"任务总数: {filtered_count}")
    if filtered_count == 0: sys.exit(0)

    header_cols = ["Protein_ID", "Res_ID", "Res_Name", "SASA", "RSA", 
                   "Charge_PQR", "Charge_Avg_Theory", "Meta_Min_PH", "Meta_Max_PH"]
    if COL_ORGANISM: header_cols.append("Meta_Organism")
    
    write_header = not os.path.exists(OUTPUT_CSV)
    
    with open(OUTPUT_CSV, "a" if not write_header else "w") as f_out:
        if write_header: f_out.write(",".join(header_cols) + "\n")
            
        processed_count = 0
        for idx, row in df.iterrows():
            processed_count += 1
            entry_id = str(row[COL_ENTRY]).strip()
            
            meta_min = row.get('temp_min_ph', '')
            meta_max = row.get('temp_max_ph', '')
            meta_org = str(row[COL_ORGANISM]).replace(",", ";") if (COL_ORGANISM and COL_ORGANISM in df.columns) else "N/A"
            if pd.isna(meta_org): meta_org = "N/A"

            print(f"[{processed_count}/{filtered_count}] 处理 {entry_id} ...", end="", flush=True)

            tmp_pdb = None
            tmp_pqr = None
            
            try:
                # 1. 下载
                tmp_pdb, final_id = download_alphafold_structure(entry_id)
                if not tmp_pdb:
                    print(f" [跳过] 无结构")
                    logging.info(f"{entry_id}: 未找到结构")
                    record_failed_id(entry_id) # 记录失败
                    time.sleep(SLEEP_ON_FAIL)
                    continue

                # 2. PDB2PQR
                tmp_pqr = run_pdb2pqr_temp(tmp_pdb, TARGET_PH)
                if not tmp_pqr:
                    print(f" [错误] PDB2PQR 失败")
                    logging.error(f"{entry_id}: PDB2PQR 计算失败")
                    record_failed_id(entry_id) # 记录失败
                    time.sleep(SLEEP_ON_FAIL)
                    continue

                # 3. 分析
                res_data = analyze_structure(tmp_pdb, tmp_pqr, TARGET_PH)
                if res_data:
                    lines = []
                    for d in res_data:
                        line_items = [
                            entry_id, str(d['res_id']), d['res_name'],
                            f"{d['sasa']:.2f}", f"{d['rsa']:.3f}",
                            f"{d['charge_pqr']:.3f}", f"{d['charge_avg']:.3f}",
                            str(meta_min), str(meta_max)
                        ]
                        if COL_ORGANISM: line_items.append(meta_org)
                        lines.append(",".join(line_items))
                    f_out.write("\n".join(lines) + "\n")
                    f_out.flush()
                    print(f" 完成")
                else:
                    print(f" [空结果]")
                    logging.warning(f"{entry_id}: 分析结果为空")
                    record_failed_id(entry_id) # 记录失败

            except Exception as e:
                print(f" [异常]")
                logging.error(f"{entry_id}: 异常: {e}")
                record_failed_id(entry_id) # 记录失败
                time.sleep(SLEEP_ON_FAIL)
                
            finally:
                if tmp_pdb and os.path.exists(tmp_pdb): os.remove(tmp_pdb)
                if tmp_pqr and os.path.exists(tmp_pqr): os.remove(tmp_pqr)

    print("\n完成！")