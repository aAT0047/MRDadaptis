import os
import subprocess

def run_breakdancer(bam_file, output_dir):
    # Step 1: 生成 BreakDancer 配置文件
    cfg_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}.cfg")
    bam2cfg_cmd = f"bam2cfg.pl -g -h -q 5 {bam_file} > {cfg_file}"
    subprocess.run(bam2cfg_cmd, shell=True)
    base_name = os.path.basename(bam_file).replace('.bam', '').split('_')[0]
    histogram_base_name = f"{os.path.basename(bam_file)}.{base_name}.insertsize_histogram"
    histogram_txt = os.path.join("/data/home/std_12/SVfolder", histogram_base_name)
    histogram_png = f"{histogram_txt}.png"  # 直方图 PNG 文件带有 .png 扩展名
    os.remove(histogram_png)  # 删除插入片段直方图 PNG 文件
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
    output_vcf = os.path.join(output_dir, f"{os.path.basename(bam_file)}.breakdancer_max.vcf")
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



bam_file = "/data/home/std_12/ShiHe/IGbam/PA235W0733-K1CFE70XNF1-L000KBKF_BCL2_1.bam"
output_dir = "/data/home/std_12"

vcf_file = run_breakdancer(bam_file, output_dir)
if vcf_file:
    print(f"VCF file created: {vcf_file}")
else:
    print("No CTX or ITX variations were found.")
