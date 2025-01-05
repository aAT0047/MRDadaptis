import os
import subprocess
import re
import shutil
from multiprocessing import Pool

from config import config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam

# 定义一个函数用于运行指定路径下的任务
def process_bam_files_for_config(config):
    input_dir = config["input_dir"]
    output_base_dir = config["output_base_dir"]

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
    tools = ["delly", "pindel", "breakdancer", "svaba"]
    tool_dirs = {tool: f"{output_base_dir}/{tool}" for tool in tools}
    for dir in tool_dirs.values():
        os.makedirs(dir, exist_ok=True)
    return tool_dirs


# 获取所有BAM文件
def get_bam_files(input_dir):
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    return bam_files



def mismatch(bam_file, output_dir):
        # 定义中间文件的路径
        
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 定义文件路径
    bam_chromosomes = os.path.join(output_dir, "bam_chromosomes.txt")
    ref_chromosomes = os.path.join(output_dir, "ref_chromosomes.txt")
    mismatched_chromosomes = os.path.join(output_dir, "mismatched_chromosomes.txt")

    # 步骤1：提取BAM文件中的染色体列表
    print("提取BAM文件中的染色体列表...")
    with open(bam_chromosomes, 'w') as bam_chr_file:
        idxstats_cmd = f"samtools idxstats {bam_file} | cut -f1"
        subprocess.run(idxstats_cmd, shell=True, stdout=bam_chr_file, check=True)

    # 步骤2：提取参考基因组中的染色体列表
    print("提取参考基因组中的染色体列表...")
    # 确保参考基因组已建立索引
    if not os.path.exists(f"{reference_genome}.fai"):
        print("参考基因组索引不存在，正在创建索引...")
        faidx_cmd = f"samtools faidx {reference_genome}"
        subprocess.run(faidx_cmd, shell=True, check=True)

    # 提取参考基因组的染色体名称
    with open(ref_chromosomes, 'w') as ref_chr_file:
        cut_cmd = f"cut -f1 {reference_genome}.fai"
        subprocess.run(cut_cmd, shell=True, stdout=ref_chr_file, check=True)

    # 步骤3：比较染色体列表，找出不匹配的染色体
    print("比较染色体列表，找出不匹配的染色体...")
    comm_cmd = f"comm -3 <(sort {bam_chromosomes}) <(sort {ref_chromosomes}) > {mismatched_chromosomes}"
    
    # 使用 bash 执行 comm 命令进行排序并比较
    subprocess.run(f"bash -c '{comm_cmd}'", shell=True, check=True)

    print(f"不匹配的染色体已保存至: {mismatched_chromosomes}")

def run_delly(bam_file, output_dir):
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
        # 运行Delly检测染色体易位
        print("运行Delly检测染色体易位（TRA）...")
        output_vcf = os.path.join(output_dir, f"{bam_prefix}.vcf")
        delly_cmd = ["delly", "call", "-t", "TRA", "-g", reference_genome, "-o", output_vcf, "-f", "bcf", filtered_bam]
        subprocess.run(delly_cmd, check=True)

        # 转换BCF文件为VCF文件

        # 检查是否存在 BCF 文件，存在则进行转换
        if os.path.exists(output_vcf):
            vcf_converted = os.path.join(output_dir, f"{bam_prefix}.vcf")
            bcftools_cmd = ["bcftools", "view", output_vcf, "-o", vcf_converted, "-O", "v"]
            subprocess.run(bcftools_cmd, check=True)
            print(f"已转换为可读VCF文件，路径：{vcf_converted}")
        else:
            print(f"{output_vcf} 文件不存在，跳过转换步骤。")



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

def run_pindel(bam_file, output_dir, insert_size=300):
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 提取 BAM 文件名（不包含路径和扩展名）
    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]

    # 定义输出文件前缀
    output_prefix = os.path.join(output_dir, bam_basename)
    sample_label = f"{bam_basename}_sample"
    # 定义合并后的 VCF 文件
    merged_vcf = os.path.join(output_dir, f"{bam_basename}.vcf.gz")

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
    # merge_cmd = f"bcftools merge {' '.join(vcf_files)} -o {merged_vcf} "
    merge_cmd = f"bcftools merge --force-samples {' '.join(vcf_files)} -o {merged_vcf}"

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
    print('finish*************************************')
    return merged_vcf

def run_svaba(bam_file, output_dir, cores=1, reference_genome="/data/home/std_12/hs37d5.fa"):
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
    os.rename(svaba_vcf, new_svaba_vcf)
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


def run_breakdancer(bam_file, output_dir):
    # Step 1: 生成 BreakDancer 配置文件
    cfg_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}.cfg")
    bam2cfg_cmd = f"bam2cfg.pl -g -h -q 5 {bam_file} > {cfg_file}"
    subprocess.run(bam2cfg_cmd, shell=True)
    base_name = os.path.basename(bam_file).replace('.bam', '').split('_')[0]
    histogram_base_name = f"{os.path.basename(bam_file)}.{base_name}.insertsize_histogram"
    histogram_txt = os.path.join("/data/home/std_12/SVfolder", histogram_base_name)
    histogram_png = f"{histogram_txt}.png"  # 直方图 PNG 文件带有 .png 扩展名
    # 检查 PNG 文件是否存在
    if os.path.exists(histogram_png):
        os.remove(histogram_png)  # 删除插入片段直方图 PNG 文件
    # 检查 TXT 文件是否存在
    if os.path.exists(histogram_txt):
        os.remove(histogram_txt)  # 删除插入片段直方图文本文件
    # Step 2: 运行 BreakDancer 检测结构变异
    output_txt = os.path.join(output_dir, f"{os.path.basename(bam_file)}.breakdancer.out")
    breakdancer_cmd = f"breakdancer-max {cfg_file} > {output_txt}"
    subprocess.run(breakdancer_cmd, shell=True)
    
    # Step 3: 查找所有 CTX 或最大的 ITX
    max_itx_size = -1
    max_itx_variation = None
    max_ctx_size = -1
    max_ctx_variation = None

    with open(output_txt, 'r') as infile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue  # 跳过注释行和空行
            columns = line.strip().split()
            sv_type = columns[6]  # 变异类型
            size = abs(int(columns[7]))  # 变异大小
            
            if sv_type == "CTX":
                if size > max_ctx_size:
                    max_ctx_size = size
                    max_ctx_variation = columns  # 保留最大的 CTX 变异
            elif sv_type == "ITX" and size > max_itx_size:
                max_itx_size = size  # 保留最大的 ITX 变异
                max_itx_variation = columns

    # Step 4: 选择 CTX 或最大的 ITX
    selected_variation = max_ctx_variation if max_ctx_variation else max_itx_variation

    # 如果没有找到 CTX 或 ITX，返回 None
    if not selected_variation:
        print("No valid CTX or ITX variations found.")
        return None

    # Step 5: 将选择的变异输出为 VCF 格式
    output_vcf = os.path.join(output_dir, f"{os.path.basename(bam_file)}.vcf")
    with open(output_vcf, 'w') as vcf:
        # 写入 VCF 头部
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=BreakDancer\n")
        vcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        vcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n")
        vcf.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
        vcf.write("##INFO=<ID=SU,Number=1,Type=Float,Description=\"Frequency(%) of supporting reads\">\n")
        vcf.write("##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Mean Depth Dedup\">\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        chrom_start = selected_variation[0]   # 染色体1
        pos_start = selected_variation[1]     # 起始位置1
        chrom_end = selected_variation[3]     # 染色体2
        pos_end = selected_variation[4]       # 结束位置2
        sv_type = selected_variation[6]       # 变异类型
        sv_size = selected_variation[7]       # 变异大小
        support_reads = selected_variation[9]  # 支持的读段数
        unique_id = "CTX_Event"  # 可以生成唯一事件ID

        if sv_type == "CTX":
            # 根据你提供的格式生成 VCF 的 ALT 和 INFO 字段
            ref = 'N'
            alt_1 = f'[{chrom_end}:{pos_end}[N'
            alt_2 = f'N[{chrom_start}:{pos_start}]'
            qual = '.'
            filter_ = 'PASS'
            info_1 = (f'SVTYPE={sv_type};STRANDS=--:503;EVENT={unique_id};MATEID=1;'
                      f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={support_reads};'
                      f'PE=0;SR=20;END={pos_end}')
            info_2 = (f'SVTYPE={sv_type};STRANDS=--:503;EVENT={unique_id};MATEID=2;'
                      f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={support_reads};'
                      f'PE=0;SR=20;END={pos_start};SECONDARY')
            format_ = 'GT:SU:PE:SR'
            sample_data = f'./.:{support_reads}:0:20'

            # 写入 VCF 内容
            vcf.write(f"{chrom_start}\t{pos_start}\t.\t{ref}\t{alt_1}\t{qual}\t{filter_}\t{info_1}\t{format_}\t{sample_data}\n")
            vcf.write(f"{chrom_end}\t{pos_end}\t.\t{ref}\t{alt_2}\t{qual}\t{filter_}\t{info_2}\t{format_}\t{sample_data}\n")

        elif sv_type == "ITX":
            # ITX 格式处理
            ref = 'N'
            alt = f"<{sv_type}>"
            qual = '.'
            filter_ = 'PASS'
            info_field = f"SVTYPE={sv_type};END={pos_end};SVLEN={sv_size};SUPPORT={support_reads}"
            format_ = 'GT:SU:PE:SR'
            sample_data = f'./.:{support_reads}:0:20'

            # 写入 VCF 内容
            vcf.write(f"{chrom_start}\t{pos_start}\t.\t{ref}\t{alt}\t{qual}\t{filter_}\t{info_field}\t{format_}\t{sample_data}\n")
    
    # Step 6: 删除中间文件，保留最终的 VCF 文件
    os.remove(cfg_file)  # 删除配置文件
    os.remove(output_txt)  # 删除 BreakDancer 输出文件

    
    return output_vcf

# 定义用于处理每个BAM文件的函数
def process_bam_file(bam_file, tool_dirs):
    bam_file_path = os.path.join(input_dir, bam_file)
    # 生成BAM索引文件（如果不存在）
    generate_bam_index(bam_file_path)
    

    
    # 运行Delly并过滤VCF
    # delly_vcf = run_delly(bam_file_path, tool_dirs['delly'])
    # filter_vcf_by_evidence(delly_vcf)
    # # 运行Pindel并过滤VCF
    # pindel_vcf = run_pindel(bam_file_path, tool_dirs['pindel'])
    # filter_vcf_by_evidence(pindel_vcf)

    # 运行SvABA并过滤VCF
    svaba_vcf= run_svaba(bam_file_path, tool_dirs['svaba'])
    filter_vcf_by_evidence(svaba_vcf)
   
    # 运行Breakdancer并过滤VCF
    # breakdancer_vcf = run_breakdancer(bam_file_path, tool_dirs['breakdancer'])
    # filter_vcf_by_evidence(breakdancer_vcf)



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

    # 定义进程池的大小，根据CPU核心数量进行设置
    pool_size = 10
    
    with Pool(pool_size) as pool:
        pool.starmap(process_bam_file, [(bam_file, tool_dirs) for bam_file in bam_files])


from tqdm import tqdm

# 脚本入口
if __name__ == '__main__':
    # configs = [config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam]
    configs = [config_ALK_RET_ROS1bam]
    bam_files = get_bam_files(input_dir)
    # 使用 tqdm 显示进度条，设置 total 参数为配置集的长度
    with tqdm(total=len(configs), desc="处理配置集", unit="bam_files") as pbar:
        # 循环处理不同的配置集
        for config in configs:
            print(f"处理 {config['input_dir']} 下的 BAM 文件")
            process_bam_files_for_config(config)
            
            # 每次处理完一个 config，更新进度条
            pbar.update(1)

    print("所有路径处理完成。")
