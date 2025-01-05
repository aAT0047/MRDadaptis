import os
import subprocess

def run_pindel(bam_file, reference_genome, output_dir="/data/home/std_12/ShiHe/results/sarcomabam", insert_size=300, sample_label="sample1"):
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 提取 BAM 文件名（不包含路径和扩展名）
    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]

    # 定义输出文件前缀
    output_prefix = os.path.join(output_dir, bam_basename)

    # 定义合并后的 VCF 文件
    merged_vcf = os.path.join(output_dir, f"{bam_basename}_merged.vcf.gz")

    # 创建 Pindel 配置文件
    config_file = os.path.join(output_dir, f"{bam_basename}_pindel_config.txt")
    with open(config_file, 'w') as config:
        config.write(f"{bam_file} {insert_size} {sample_label}\n")

    # 运行 Pindel
    pindel_cmd = f"pindel -f '{reference_genome}' -i {config_file} -o {output_prefix}"
    subprocess.run(pindel_cmd, shell=True, check=True)

    # 处理的后缀列表
    suffixes = ['_D', '_INV', '_TD', '_BP', '_SI', '_LI', '_RP', '_CloseEndMapped', '_INT_final']
    vcf_files = []

    # 转换所有结果文件为 VCF，压缩并索引
    for suffix in suffixes:
        result_file = f"{output_prefix}{suffix}"
        vcf_file = f"{output_prefix}{suffix}.vcf"
        vcf_gz_file = f"{vcf_file}.gz"
        
        # 转换为 VCF
        vcf_cmd = f"pindel2vcf -p {result_file} -r {reference_genome} -R human -d 20240101 -v {vcf_file}"
        subprocess.run(vcf_cmd, shell=True, check=True)
        vcf_files.append(vcf_gz_file)

        # 压缩 VCF 文件
        bgzip_cmd = f"bgzip {vcf_file}"
        subprocess.run(bgzip_cmd, shell=True, check=True)

        # 创建索引
        index_cmd = f"bcftools index {vcf_gz_file}"
        subprocess.run(index_cmd, shell=True, check=True)

        print(f"已生成压缩并索引的 VCF 文件：{vcf_gz_file}")

    # 合并多个 VCF.gz 文件为一个
    merge_cmd = f"bcftools merge {' '.join(vcf_files)} -o {merged_vcf} "
    subprocess.run(merge_cmd, shell=True, check=True)
    print(f"所有 VCF 文件已合并为：{merged_vcf}")

    # 删除中间文件
    intermediate_files = [config_file] + [f"{output_prefix}{suffix}" for suffix in suffixes]
    for file in intermediate_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"已删除中间文件: {file}")
    
    # 删除生成的单独的 VCF.gz 文件及其索引
    for vcf_file in vcf_files:
        if os.path.exists(vcf_file):
            os.remove(vcf_file)
            print(f"已删除单独的 VCF.gz 文件: {vcf_file}")
        
        index_file = f"{vcf_file}.csi"
        if os.path.exists(index_file):
            os.remove(index_file)
            print(f"已删除索引文件: {index_file}")

    return merged_vcf


bam_file = "/data/home/std_12/ShiHe/sarcomabam/TB233V0923-D2ABS30XNF1-L000_STAT6_1.bam"
reference_genome = "/data/home/std_12/GRCh38/GRCh38.fa"

# 您可以显式传递 output_dir 或使用默认路径
run_pindel(bam_file, reference_genome)
