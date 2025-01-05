import subprocess
import os

def run_delly_with_filtering(bam_file, reference_genome, output_dir):
    """
    提取BAM文件和参考基因组中的染色体列表，过滤BAM文件中不匹配的染色体，运行Delly，并清理临时文件。

    参数：
    bam_file (str): 输入的BAM文件路径。
    reference_genome (str): 参考基因组的FASTA文件路径。
    output_dir (str): 输出文件保存的目录。

    返回：
    str: Delly生成的VCF文件路径。
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 提取BAM文件的基本信息
    bam_basename = os.path.basename(bam_file)
    bam_prefix = os.path.splitext(bam_basename)[0]

    # 定义中间文件的路径
    bam_chromosomes = os.path.join(output_dir, f"{bam_prefix}_bam_chromosomes.txt")
    ref_chromosomes = os.path.join(output_dir, f"{bam_prefix}_ref_chromosomes.txt")
    mismatched_chromosomes = os.path.join(output_dir, f"{bam_prefix}_mismatched_chromosomes.txt")
    filtered_bam = os.path.join(output_dir, f"{bam_prefix}_filtered.bam")

    try:
        # 步骤1：提取BAM文件中的染色体列表
        print(f"提取 {bam_file} 中的染色体列表...")
        idxstats_cmd = f"samtools idxstats {bam_file} | cut -f1"
        with open(bam_chromosomes, 'w') as bam_chr_file:
            subprocess.run(idxstats_cmd, shell=True, stdout=bam_chr_file, check=True)

        # 步骤2：提取参考基因组中的染色体列表
        print("提取参考基因组中的染色体列表...")
        if not os.path.exists(f"{reference_genome}.fai"):
            print(f"参考基因组索引不存在，正在创建索引...")
            subprocess.run(f"samtools faidx {reference_genome}", shell=True, check=True)
        cut_cmd = f"cut -f1 {reference_genome}.fai"
        with open(ref_chromosomes, 'w') as ref_chr_file:
            subprocess.run(cut_cmd, shell=True, stdout=ref_chr_file, check=True)

        # 步骤3：比较染色体列表，找出不匹配的染色体
        print("比较染色体列表，找出不匹配的染色体...")
        comm_cmd = f"comm -3 <(sort {bam_chromosomes}) <(sort {ref_chromosomes}) > {mismatched_chromosomes}"
        subprocess.run(f"bash -c '{comm_cmd}'", shell=True, check=True)

        # 步骤4：过滤BAM文件，排除不匹配的染色体
        print("开始过滤BAM文件，排除不匹配的染色体...")
        filter_cmd = f"samtools view -h {bam_file} | grep -v -F -f {mismatched_chromosomes} | samtools view -b -o {filtered_bam}"
        subprocess.run(f"bash -c '{filter_cmd}'", shell=True, check=True)

        # 步骤5：为过滤后的BAM文件创建索引
        print("正在为过滤后的BAM文件建立索引...")
        subprocess.run(f"samtools index {filtered_bam}", shell=True, check=True)

        # 步骤6：运行Delly
        print("运行Delly...")
        output_vcf = os.path.join(output_dir, f"{bam_prefix}.delly.vcf")
        delly_cmd = ["delly", "call","-t", "TRA", "-g", reference_genome, "-o", output_vcf, filtered_bam]
        subprocess.run(delly_cmd, check=True)
        print(f"Delly分析已完成，生成的VCF文件位于：{output_vcf}")

        # 步骤7：清理临时文件
        print("清理临时文件...")
        temp_files = [bam_chromosomes, ref_chromosomes, mismatched_chromosomes]
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"已删除临时文件: {temp_file}")

        # 删除过滤后的BAM文件及其索引文件
        if os.path.exists(filtered_bam):
            os.remove(filtered_bam)
            print(f"已删除过滤后的BAM文件: {filtered_bam}")
        filtered_bam_index = filtered_bam + '.bai'
        if os.path.exists(filtered_bam_index):
            os.remove(filtered_bam_index)
            print(f"已删除过滤后的BAM索引文件: {filtered_bam_index}")

    except subprocess.CalledProcessError as e:
        print(f"操作时出错: {e}")
        return None

    return output_vcf

# 使用示例
bam_file = "/data/home/std_12/ShiHe/sarcomabam/FN24270111-L1ABAW2GXF1-L000KA46_EWSR1_1.bam"
reference_genome = "/data/home/std_12/GRCh38/GRCh38.fa"
output_dir = "/data/home/std_12/output"

vcf_file = run_delly_with_filtering(bam_file, reference_genome, output_dir)
print(f"最终生成的VCF文件：{vcf_file}")
