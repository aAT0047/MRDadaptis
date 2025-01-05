import os
import subprocess
import csv
from multiprocessing import Pool
from evaf1 import *
# ------------------------------
# 全局变量定义和路径设置
# ------------------------------
base_dir = "/data/home/std_12/ShiHe"
input_dir = f"{base_dir}/sarcomabam"
output_base_dir = f"{base_dir}/results"

# 用于比对的真值 VCF 文件目录
vcf_dir = f"{base_dir}/sarcomavcf"

# 定义要运行的工具列表（不包括 Manta）
tools = ["delly", "pindel", "breakdancer", "svaba"]
tool_dirs = {tool: f"{output_base_dir}/{tool}" for tool in tools}
for dir in tool_dirs.values():
    os.makedirs(dir, exist_ok=True)

# 创建用于存储评估结果的目录
evaluation_dir = f"{output_base_dir}/evaluation"
os.makedirs(evaluation_dir, exist_ok=True)

# ------------------------------
# 获取所有 BAM 文件
# ------------------------------
def get_bam_files(input_dir):
    """获取输入目录中所有的 BAM 文件列表。"""
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    return bam_files

bam_files = get_bam_files(input_dir)

# ------------------------------
# 评估函数定义
# ------------------------------
def filter_vcf_by_evidence(vcf_file):
    """
    根据证据过滤 VCF 文件。
    这里需要根据实际需求实现具体的过滤逻辑。
    示例：仅保留 PASS 过滤的变异。
    """
    filtered_vcf = vcf_file.replace('.vcf', '.filtered.vcf')
    with open(vcf_file, 'r') as infile, open(filtered_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                if 'PASS' in fields[6]:
                    outfile.write(line)
    return filtered_vcf

def evamain(ground_truth_vcf, detected_vcf):
    """
    评估检测的 VCF 文件与真值 VCF 文件的匹配情况。
    
    参数：
    - ground_truth_vcf: 真值 VCF 文件路径
    - detected_vcf: 检测到的 VCF 文件路径
    
    返回：
    - TP, FP, FN, F1_score, Precision, Recall
    """
    # 简单示例：基于变异数量进行评估。
    # 实际应用中应使用更复杂的匹配逻辑，如基于变异位置和类型的匹配。
    
    def parse_vcf(vcf):
        variants = set()
        with open(vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = fields[1]
                var_type = fields[2]  # 假设第三列表示变异类型
                variants.add((chrom, pos, var_type))
        return variants

    truth_variants = parse_vcf(ground_truth_vcf)
    detected_variants = parse_vcf(detected_vcf)
    
    TP = len(truth_variants & detected_variants)
    FP = len(detected_variants - truth_variants)
    FN = len(truth_variants - detected_variants)
    
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return TP, FP, FN, f1_score, precision, recall

# ------------------------------
# 运行工具的函数定义
# ------------------------------
def run_tool(tool, bam_file, output_dir):
    """
    运行指定的工具。
    
    参数：
    - tool: 工具名称（delly, pindel, breakdancer, svaba）
    - bam_file: BAM 文件路径
    - output_dir: 工具的输出目录
    
    返回：
    - output_vcf: 生成的 VCF 文件路径
    """
    sample_name = os.path.splitext(os.path.basename(bam_file))[0]
    output_vcf = os.path.join(output_dir, f"{sample_name}.{tool}.vcf")
    
    if tool == "delly":
        delly_cmd = f"delly call -o {output_vcf} {bam_file}"
        subprocess.run(delly_cmd, shell=True, check=True)
    elif tool == "pindel":
        # Pindel 的输出通常不是 VCF 格式，需要转换，这里假设已转换为 VCF
        pindel_cmd = f"pindel -f {reference_genome} -i {bam_file} -o {output_vcf}"
        subprocess.run(pindel_cmd, shell=True, check=True)
    elif tool == "breakdancer":
        # Breakdancer 生成的是 DEL 文件，需转换为 VCF，这里假设已转换为 VCF
        breakdancer_cmd = f"breakdancer -o {output_vcf} {bam_file}"
        subprocess.run(breakdancer_cmd, shell=True, check=True)
    elif tool == "svaba":
        svaba_cmd = f"svaba run -t {bam_file} -a {output_vcf}"
        subprocess.run(svaba_cmd, shell=True, check=True)
    else:
        print(f"未知的工具: {tool}")
        return None
    
    return output_vcf

# ------------------------------
# 评估并记录结果的函数
# ------------------------------
def evaluate_and_record(tool, sample_name, detected_vcf, csv_writers):
    """
    对检测到的 VCF 进行评估，并将结果记录到 CSV 文件中。
    
    参数：
    - tool: 工具名称
    - sample_name: 样本名称
    - detected_vcf: 检测到的 VCF 文件路径
    - csv_writers: 字典，键为工具名称，值为对应的 CSV 写入器
    """
    # 获取对应的真值 VCF 文件路径
    ground_truth_vcf = os.path.join(vcf_dir, f"{sample_name}.vcf")
    
    if not os.path.exists(ground_truth_vcf):
        print(f"真值 VCF 文件不存在: {ground_truth_vcf}")
        return
    
    # 过滤 VCF 文件
    filtered_vcf = filter_vcf_by_evidence(detected_vcf)
    
    # 评估
    TP, FP, FN, f1_score, precision, recall = evamain(ground_truth_vcf, filtered_vcf)
    
    # 写入 CSV
    csv_writers[tool].writerow({
        "Sample": sample_name,
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "Precision": precision,
        "Recall": recall,
        "F1_Score": f1_score
    })

# ------------------------------
# 处理每个 BAM 文件的函数
# ------------------------------
def process_bam_file(bam_file):
    """
    处理单个 BAM 文件，运行所有工具并进行评估。
    
    参数：
    - bam_file: BAM 文件名
    """
    sample_name = os.path.splitext(os.path.basename(bam_file))[0]
    bam_file_path = os.path.join(input_dir, bam_file)
    
    results = []
    
    # 运行所有工具
    for tool in tools:
        output_dir = tool_dirs[tool]
        try:
            output_vcf = run_tool(tool, bam_file_path, output_dir)
            if output_vcf and os.path.exists(output_vcf):
                results.append((tool, sample_name, output_vcf))
            else:
                print(f"{tool} 未生成有效的 VCF 文件: {output_vcf}")
        except subprocess.CalledProcessError as e:
            print(f"运行 {tool} 时出错: {e}")
    
    # 评估并记录结果
    for tool, sample_name, detected_vcf in results:
        evaluate_and_record(tool, sample_name, detected_vcf, csv_writers)

# ------------------------------
# 主执行块，使用多核并行处理
# ------------------------------
if __name__ == '__main__':
    # 定义进程池的大小，基于可用的 CPU 核心数量
    pool_size = os.cpu_count() or 1  # 如果无法获取，则默认为 1
    
    print(f"开始使用 {pool_size} 个并行进程处理 BAM 文件和评估结果...")
    
    # 初始化 CSV 写入器
    csv_files = {}
    csv_writers = {}
    for tool in tools:
        csv_path = os.path.join(evaluation_dir, f"{tool}_evaluation.csv")
        csv_file = open(csv_path, 'w', newline='')
        csv_writer = csv.DictWriter(csv_file, fieldnames=["Sample", "TP", "FP", "FN", "Precision", "Recall", "F1_Score"])
        csv_writer.writeheader()
        csv_files[tool] = csv_file
        csv_writers[tool] = csv_writer
    
    try:
        with Pool(pool_size) as pool:
            # 并行处理所有 BAM 文件
            pool.map(process_bam_file, bam_files)
    finally:
        # 关闭所有 CSV 文件
        for csv_file in csv_files.values():
            csv_file.close()
    
    # 合并所有 CSV 文件到一个汇总的 CSV
    summary_csv_path = os.path.join(evaluation_dir, "summary_evaluation.csv")
    with open(summary_csv_path, 'w', newline='') as summary_csv:
        fieldnames = ["Tool", "Sample", "TP", "FP", "FN", "Precision", "Recall", "F1_Score"]
        writer = csv.DictWriter(summary_csv, fieldnames=fieldnames)
        writer.writeheader()
        
        for tool in tools:
            tool_csv_path = os.path.join(evaluation_dir, f"{tool}_evaluation.csv")
            if not os.path.exists(tool_csv_path):
                print(f"工具 {tool} 的评估 CSV 文件不存在: {tool_csv_path}")
                continue
            with open(tool_csv_path, 'r') as tool_csv:
                reader = csv.DictReader(tool_csv)
                for row in reader:
                    summary_row = {"Tool": tool}
                    summary_row.update(row)
                    writer.writerow(summary_row)
    
    print("所有评估结果已完成并合并到 summary_evaluation.csv。")
