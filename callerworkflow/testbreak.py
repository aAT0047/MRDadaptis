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



def run_breakdancer(bam_file, output_dir):
    # Step 1: 生成 BreakDancer 配置文件
    cfg_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}.cfg")
    bam2cfg_cmd = f"bam2cfg.pl -g -h -q 5 {bam_file} > {cfg_file}"
    subprocess.run(bam2cfg_cmd, shell=True)

    # 获取基本文件名信息
    base_name = os.path.basename(bam_file).replace('.bam', '').split('_')[0]
    histogram_base_name = f"{os.path.basename(bam_file)}.{base_name}.insertsize_histogram"
    histogram_txt = os.path.join("/data/home/std_12/SVfolder", histogram_base_name)
    histogram_png = f"{histogram_txt}.png"  # 直方图 PNG 文件

    # 检查并删除 PNG 和 TXT 文件
    if os.path.exists(histogram_png):
        os.remove(histogram_png)
    if os.path.exists(histogram_txt):
        os.remove(histogram_txt)
    if os.path.exists(os.path.join("/data/home/std_12", f"{histogram_base_name}.png")):
        os.remove(os.path.join("/data/home/std_12", f"{histogram_base_name}.png"))
    if os.path.exists(os.path.join("/data/home/std_12", histogram_base_name)):
        os.remove(os.path.join("/data/home/std_12", histogram_base_name))

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

            # 检查列数，避免超出索引
            if len(columns) < 8:
                print(f"Skipping malformed line: {line.strip()}")
                continue

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
        unique_id = "CTX_Event"  # 生成唯一事件ID

        if sv_type == "CTX":
            # 生成 VCF 格式的 ALT 和 INFO 字段
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

    # # 运行SvABA并过滤VCF
    # svaba_vcf= run_svaba(bam_file_path, tool_dirs['svaba'])
    # filter_vcf_by_evidence(svaba_vcf)
   
    # 运行Breakdancer并过滤VCF
    breakdancer_vcf = run_breakdancer(bam_file_path, tool_dirs['breakdancer'])
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
    pool_size = 20
    
    with Pool(pool_size) as pool:
        pool.starmap(process_bam_file, [(bam_file, tool_dirs) for bam_file in bam_files])

# 脚本入口
if __name__ == '__main__':
    configs = [config_IGbam, config_ALK_RET_ROS1bam, config_sarcomabam]
    
    # 循环处理不同的配置集
    for config in configs:
        print(f"处理 {config['input_dir']} 下的 BAM 文件")
        process_bam_files_for_config(config)

    print("所有路径处理完成。")