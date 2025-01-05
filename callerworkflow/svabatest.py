import os
import subprocess
import vcfpy
import re
import shutil
from config import config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam

# 定义生成BAM索引的函数
def generate_bam_index(bam_file):
    bam_index = bam_file + ".bai"
    if not os.path.exists(bam_index):
        print(f"索引文件 {bam_index} 不存在，正在生成...")
        subprocess.run(f"samtools index {bam_file}", shell=True)
    else:
        print(f"索引文件 {bam_index} 已存在，跳过生成。")

# 定义工具目录创建
def create_tool_dirs(output_base_dir):
    tools = ["svaba"]
    tool_dirs = {tool: f"{output_base_dir}/{tool}" for tool in tools}
    for dir in tool_dirs.values():
        os.makedirs(dir, exist_ok=True)
    return tool_dirs

# 获取所有BAM文件
def get_bam_files(input_dir):
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    return bam_files
    
def run_svaba(bam_file, output_dir, cores=10, reference_genome="/data/home/std_12/hs37d5.fa"):
    """
    运行 SvABA，生成 VCF 文件，并过滤结果，删除中间文件
    """
    # Step 1: 设置输出前缀
    output_prefix = os.path.join(output_dir, f"{os.path.basename(bam_file).replace('.bam', '')}")
    
    # Step 2: 构造 SvABA 命令
    svaba_cmd = f"svaba run -t {bam_file} -p {cores} -a {output_prefix} -G {reference_genome}"
    print(f"Running command: {svaba_cmd}")
    subprocess.run(svaba_cmd, shell=True)
    
    # Step 3: 获取 SvABA 生成的 VCF 文件
    svaba_vcf = f"{output_prefix}.svaba.sv.vcf"
    # 定义新的文件名
    new_svaba_vcf = f"{output_prefix}.vcf"
   # os.rename(svaba_vcf, new_svaba_vcf)
        # 检查 VCF 文件是否存在，存在则重命名
    if os.path.exists(svaba_vcf):
        os.rename(svaba_vcf, new_svaba_vcf)
        print(f"Renamed {svaba_vcf} to {new_svaba_vcf}")
    else:
        print(f"File {svaba_vcf} not found, skipping renaming...")
        # 删除不需要的文件
    files_to_remove = [
        f"{output_prefix}.alignments.txt.gz",
        f"{output_prefix}.bps.txt.gz",
        f"{output_prefix}.contigs.bam",
        f"{output_prefix}.discordant.txt.gz",
        f"{output_prefix}.log",
        f"{output_prefix}.svaba.indel.vcf",
        f"{output_prefix}.svaba.unfiltered.indel.vcf",
        f"{output_prefix}.svaba.unfiltered.sv.vcf"
    ]

    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")
        else:
            print(f"File {file} not found, skipping...")

    return new_svaba_vcf

def filter_vcf_by_evidence(vcf_file):
    """
    过滤 VCF 文件，只保留支持证据最多的变异对，并返回这些变异行。
    
    :param vcf_file: 输入 VCF 文件路径
    :return: 支持证据最多的变异行列表
    """
    temp_vcf_file = vcf_file + ".tmp"

    # 读取 VCF 文件
    with open(vcf_file, 'r') as f:
        lines = f.readlines()

    # 提取变异信息
    variants = {}
    header = []
    for line in lines:
        if line.startswith("#"):
            header.append(line)
            continue
        parts = line.strip().split("\t")
        info = parts[7]
                # 提取 MAPQ, NM, ID 和 MATEID
        mapq_match = re.search(r'MAPQ=(\d+)', info)
        nm_match = re.search(r'NM=(\d+)', info)
        id_match =parts[2] # 提取变异的唯一ID
        mate_id_match = re.search(r'MATEID=([^;]+)', info)  # 提取配对断点ID

        if mapq_match and nm_match and id_match and mate_id_match:
            mapq = int(mapq_match.group(1))
            nm = int(nm_match.group(1))
            var_id = id_match
            mate_id = mate_id_match.group(1)

            # 将每个变异信息存入 variants 字典，键为 ID，值为包含所有信息的元组
            variants[var_id] = (line, mapq - nm, mate_id)

    # 存储已经处理的变异对，避免重复
    processed_pairs = set()

    # 选择 MAPQ 最高且 NM 最低的变异对
    best_variant_pair = []
    highest_score = float('-inf')

    for var_id, (line, score, mate_id) in variants.items():
        # 确保变异对存在并且没有重复处理
        if mate_id in variants and (var_id, mate_id) not in processed_pairs and (mate_id, var_id) not in processed_pairs:
            # 计算该变异对的总得分
            total_score = score + variants[mate_id][1]

            # 如果当前变异对得分高于之前的最佳得分，则更新
            if total_score > highest_score:
                highest_score = total_score
                best_variant_pair = [line, variants[mate_id][0]]

            # 标记这对变异为已处理
            processed_pairs.add((var_id, mate_id))

    # 生成过滤后的临时 VCF 文件
    with open(temp_vcf_file, 'w') as f:
        for line in header:
            f.write(line)
        for variant in best_variant_pair:
            f.write(variant)

    # 用临时 VCF 文件覆盖原始 VCF 文件
    shutil.move(temp_vcf_file, vcf_file)
    # os.remove(temp_vcf_file)
    print(vcf_file)

# 主函数，处理不同的配置
def process_bam_files_for_config(config):
    global input_dir, output_base_dir, reference_genome
    input_dir = config["input_dir"]
    output_base_dir = config["output_base_dir"]
    reference_genome = '/data/home/std_12/hs37d5.fa'
    
    # 获取BAM文件列表
    bam_files = get_bam_files(input_dir)
    
    # 创建工具输出文件夹
    tool_dirs = create_tool_dirs(output_base_dir)

    # 顺序处理每个BAM文件
    for bam_file in bam_files:
        bam_file_path = os.path.join(input_dir, bam_file)
        
        # 生成BAM索引文件（如果不存在）
        generate_bam_index(bam_file_path)
        
        # 运行SvABA并过滤VCF
        svaba_vcf = run_svaba(bam_file_path, tool_dirs['svaba'])
        filter_vcf_by_evidence(svaba_vcf)
        print('finish 88888888888' ,svaba_vcf,'finish 88888888888 ,')

# 脚本入口
if __name__ == '__main__':
#    configs = [config_ALK_RET_ROS1bam,config_IGbam]
  #  configs = [config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam] 
    configs = [config_sarcomabam]
    # 循环处理不同的配置集
    for config in configs:
        print(f"处理 {config['input_dir']} 下的 BAM 文件")
        process_bam_files_for_config(config)

    print("所有路径处理完成。")
