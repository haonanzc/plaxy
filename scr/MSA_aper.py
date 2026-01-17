import os
import time
from datetime import datetime
from openpyxl import Workbook

def main():
    time_start = time.time()

    # --- 1. 获取输入文件 ---
    file_path = input("请输入fasta文件路径: ").strip()
    # 去除可能存在的引号
    file_path = file_path.strip('"').strip("'")

    if not os.path.exists(file_path):
        print(f"错误: 文件 '{file_path}' 不存在!")
        return

    # --- 2. 健壮的 FASTA 读取逻辑 ---
    # 这种方法不依赖固定的行数，能处理标准格式和多行折叠格式的FASTA
    names = []
    sequences = []
    current_seq_parts = []
    
    print("正在读取并解析 FASTA 文件...")
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line: continue # 跳过空行
            
            if line.startswith('>'):
                # 如果当前有正在累积的序列，先存起来
                if current_seq_parts:
                    sequences.append("".join(current_seq_parts))
                    current_seq_parts = []
                names.append(line[1:]) # 去掉 '>'
            else:
                current_seq_parts.append(line)
        
        # 循环结束后，保存最后一条序列
        if current_seq_parts:
            sequences.append("".join(current_seq_parts))

    # 简单校验
    if not sequences:
        print("错误: 未找到序列信息。")
        return
    
    seq_len = len(sequences[0])
    num_seqs = len(sequences)
    print(f"共读取到 {num_seqs} 条序列，序列长度为: {seq_len}")
    
    # 检查MSA对齐情况 (可选警告)
    if any(len(s) != seq_len for s in sequences):
        print("警告: 检测到序列长度不一致！这可能不是一个对齐过的(MSA)文件，统计结果可能无意义。")

    # --- 3. 初始化 Excel ---
    wb = Workbook()
    ws = wb.active
    
    # 定义表头
    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    headers = ['Pos', 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-', 
               'Seq1_aa', 'Top1_AA', 'Top1_%', 'Top2_AA', 'Top2_%', 'Top3_AA', 'Top3_%']
    ws.append(headers)

    # --- 4. 逐位点统计 (核心修改) ---
    print("正在计算位点频率...")
    
    # 获取第一条序列作为参考显示 (即便它不是WT，显示出来也有助于定位)
    main_seq = sequences[0]

    for i in range(seq_len):
        # 1. 提取当前列的所有氨基酸
        current_column = [seq[i] for seq in sequences]
        
        # 2. 基础统计
        row_data = [i + 1] # Pos (从1开始)
        
        gap_count = current_column.count('-')
        valid_aa_count = num_seqs - gap_count # 有效氨基酸总数 (作为分母)
        
        # 临时字典用于存储当前列各氨基酸的百分比，方便后续排序
        col_stats = [] 

        # 3. 计算20种氨基酸频率
        for aa in aa_list:
            count = current_column.count(aa)
            
            if valid_aa_count > 0:
                # 按照 原脚本逻辑：频率 = 该氨基酸数 / (总序列数 - Gap数)
                # 这表示"在非Gap序列中的出现频率"
                percent = count / valid_aa_count
            else:
                # 如果这一列全是 Gap，频率记为 0
                percent = 0.0
            
            row_data.append(percent)
            col_stats.append({"aa": aa, "per": percent})

        # 添加 Gap 的频率 (相对于总序列数)
        row_data.append(gap_count / num_seqs) 
        
        # 添加 Seq1 的氨基酸 (原 wt_aa)
        row_data.append(main_seq[i])

        # 4. 排序提取 Top 1/2/3
        # 按照频率降序排列
        sorted_stats = sorted(col_stats, key=lambda x: x["per"], reverse=True)
        
        # 填入 Top 3
        for k in range(3):
            if k < len(sorted_stats):
                row_data.append(sorted_stats[k]['aa'])
                row_data.append(sorted_stats[k]['per'])
            else:
                row_data.append("")
                row_data.append(0)

        # 写入 Excel
        ws.append(row_data)

    # --- 5. 处理输出路径 ---
    # 获取用户自定义路径
    user_out_dir = input("请输入保存目录 (直接回车则保存在脚本所在目录): ").strip()
    user_out_dir = user_out_dir.strip('"').strip("'")
    
    if user_out_dir and not os.path.exists(user_out_dir):
        try:
            os.makedirs(user_out_dir)
            print(f"已创建目录: {user_out_dir}")
        except Exception as e:
            print(f"无法创建目录，将保存至默认路径。错误: {e}")
            user_out_dir = ""

    now = datetime.now()
    formatted_now = now.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = f"MSA_Stats_{formatted_now}.xlsx"
    
    if user_out_dir:
        save_path = os.path.join(user_out_dir, file_name)
    else:
        save_path = os.path.join(os.getcwd(), file_name)

    try:
        wb.save(save_path)
        print("-" * 30)
        print(f"成功! 文件已保存至:\n{os.path.abspath(save_path)}")
    except Exception as e:
        print(f"保存失败 (可能是文件被占用): {e}")

    time_end = time.time()
    print(f"耗时: {time_end - time_start:.2f} 秒")

if __name__ == "__main__":
    main()