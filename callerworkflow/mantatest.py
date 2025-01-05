# -*- coding: utf-8 -*-
from __future__ import print_function  # 引入 Python 3 的 print 函数
import os
import subprocess
import shutil
import re
import gzip


# 定义参考基因组
reference_genome = '/data/home/std_12/hs37d5.fa'

# 定义Manta运行的函数
def run_manta(bam_file, output_dir):
    """
    运行Manta以检测结构变异，处理BAM文件，并将结果保存到输出目录。
    参数:
        bam_file (str): 输入BAM文件的路径
        output_dir (str): 输出文件保存的目录
    返回:
        manta_vcf (str): Manta生成的VCF文件路径
    """
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 为每个 BAM 文件生成唯一的 Manta 运行目录
    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
    run_dir = os.path.join(output_dir, "Manta_Output_{}".format(bam_basename))

    # 检查运行目录是否已经存在，如果存在则删除
    if os.path.exists(run_dir):
        print("运行目录已存在，正在删除旧目录：{}".format(run_dir))
        shutil.rmtree(run_dir)  # 删除旧的 Manta 运行目录

    os.makedirs(run_dir)

    # Step 1: 配置Manta
    print("配置Manta运行：{}".format(bam_file))
    config_cmd = "configManta.py --tumorBam {} --referenceFasta {} --runDir {}".format(bam_file, reference_genome, run_dir)
    subprocess.check_call(config_cmd, shell=True)

    # Step 2: 运行Manta的工作流程
    print("开始运行Manta工作流：{}".format(run_dir))
    workflow_cmd = "python {}/runWorkflow.py -m local -j 20".format(run_dir)
    subprocess.check_call(workflow_cmd, shell=True)

    # Manta 生成的 VCF 文件路径
    manta_vcf = os.path.join(run_dir, "results/variants/tumorSV.vcf.gz")
    print("Manta VCF 生成：{}".format(manta_vcf))
    
    return manta_vcf

# 过滤VCF文件，只保留证据最多的一对断点重排
def filter_vcf_for_best_bnd_pair(vcf_file):
    """
    过滤VCF文件，只保留支持证据最多的一对断点重排。
    
    参数:
        vcf_file (str): 输入VCF文件的路径
    返回:
        filtered_vcf (str): 过滤后的VCF文件路径
    """
    temp_vcf_file = vcf_file + ".tmp"
    best_bnd_pair = []
    highest_bnd_depth = 0

    # 打开VCF文件（支持gz压缩的VCF）
    open_func = gzip.open if vcf_file.endswith(".gz") else open
    with open_func(vcf_file, 'rt') as f:
        lines = f.readlines()

    header = []
    bnd_variants = []

    # 提取BND变异的信息
    for line in lines:
        if line.startswith("#"):
            header.append(line)  # 保存VCF头部信息
        else:
            parts = line.strip().split("\t")
            info = parts[7]
            if "SVTYPE=BND" in info:
                # 提取支持证据的深度（BND_DEPTH 字段）
                match = re.search(r'BND_DEPTH=(\d+)', info)
                if match:
                    bnd_depth = int(match.group(1))
                    bnd_variants.append((line, bnd_depth))

    # 筛选支持证据最多的一对断点重排
    for variant in bnd_variants:
        line, bnd_depth = variant
        if bnd_depth > highest_bnd_depth:
            highest_bnd_depth = bnd_depth
            best_bnd_pair = [line]  # 保存当前最高的BND变异对

    # 将结果写入临时的VCF文件，只保留头部和筛选出的变异
    with open(temp_vcf_file, 'w') as f:
        for line in header:
            f.write(line)
        for line in best_bnd_pair:
            f.write(line)

    # 用临时文件覆盖原始VCF文件
    filtered_vcf = vcf_file.replace(".vcf.gz", ".filtered.vcf")
    shutil.move(temp_vcf_file, filtered_vcf)

    print("已生成过滤后的VCF文件：{}".format(filtered_vcf))
    return filtered_vcf

# 获取所有BAM文件
def get_bam_files(input_dir):
    """
    获取指定目录下的所有 BAM 文件。
    参数:
        input_dir (str): BAM 文件所在的输入目录
    返回:
        bam_files (list): BAM 文件路径列表
    """
    bam_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.bam')]
    return bam_files


# 定义主函数处理BAM文件
def process_bam_files_for_manta(config, final_output_base_dir):
    """
    处理指定目录下的所有 BAM 文件并运行 Manta，然后将 VCF 文件移动到目标目录。
    参数:
        config (dict): 包含输入目录和输出目录的配置
        final_output_base_dir (str): 存放最终VCF文件的根目录，根据配置创建子目录
    """
    input_dir = config["input_dir"]
    output_base_dir = config["output_base_dir"]

    # 根据配置集创建对应的子目录
    config_name = os.path.basename(input_dir)  # 使用输入目录的名称作为子目录名
    final_output_dir = os.path.join(final_output_base_dir, config_name)

    # 确保最终输出目录存在
    if not os.path.exists(final_output_dir):
        os.makedirs(final_output_dir)

    # 获取BAM文件列表
    bam_files = get_bam_files(input_dir)

    # 循环处理每个BAM文件
    for bam_file in bam_files:
        print("处理BAM文件：{}".format(bam_file))
        output_dir = os.path.join(output_base_dir, os.path.basename(bam_file).replace(".bam", ""))
        
        # Step 1: 运行Manta
        manta_vcf = run_manta(bam_file, output_dir)

        # Step 2: 过滤VCF，只保留证据最多的断点重排
        filtered_vcf = filter_vcf_for_best_bnd_pair(manta_vcf)

        # Step 3: 重命名VCF文件为BAM文件的名称
        bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
        final_vcf_name = "{}.vcf".format(bam_basename)

        # 移动最终VCF文件到目标目录
        final_vcf_path = os.path.join(final_output_dir, final_vcf_name)
        shutil.move(filtered_vcf, final_vcf_path)  # 直接移动未压缩的VCF文件
        print("已移动过滤后的VCF文件到：{}".format(final_vcf_path))

# 定义配置
config_IGbam = {
    "input_dir": "/data/home/std_12/ShiHe/IGbam",
    "output_base_dir": "/data/home/std_12/Manta_Results/IGbam"
}

config_ALK_RET_ROS1bam = {
    "input_dir": "/data/home/std_12/ShiHe/ALK-RET-ROS1bam",
    "output_base_dir": "/data/home/std_12/Manta_Results/ALK-RET-ROS1bam"
}

config_sarcomabam = {
    "input_dir": "/data/home/std_12/ShiHe/sarcomabam",
    "output_base_dir": "/data/home/std_12/Manta_Results/sarcomabam"
}

# 脚本入口
if __name__ == '__main__':    
    # 指定最终VCF文件的根目录，不同配置集将会在该目录下创建子目录
    final_output_base_dir = "/data/home/std_12/Manta_Results/All_Samples_VCF"

    configs = [config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam]
    
    # 循环处理不同的配置集
    for config in configs:
        print("处理 {} 下的 BAM 文件".format(config['input_dir']))
        process_bam_files_for_manta(config, final_output_base_dir)

    print("所有BAM文件处理完成，VCF文件已按配置集分类存放。")
