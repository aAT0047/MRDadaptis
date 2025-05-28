import os
import subprocess
import pandas as pd

# 读取CSV文件并清理列名
def read_csv(file_path, encoding='latin1'):
    df = pd.read_csv(file_path, encoding=encoding)
    df.columns = [
        'Sample_ID', 'Mutation_Type', 'Mutation_Base',
        'NM_Code', 'Mutation_Form/Product', 'Chr_Start', 
        'Chr_End', 'Frequency(%)', 'Mean_Depth_Dedup', 'Is_Panel_Sample'
    ]
    return df

# 确保目录存在，不存在则创建
def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# 使用samtools裁剪BAM文件
def clip_bam(sample_id, chrom, start, end, bam_dir, output_bam):
    bam_file = os.path.join(bam_dir, f'{sample_id}.bam')  # 假设bam文件命名为Sample_ID.bam
    if not os.path.exists(bam_file):
        print(f"BAM file not found: {bam_file}")
        return
    
    # 使用samtools view裁剪BAM，并启用多线程
    region = f"{chrom}:{start}-{end}"
    cmd = ['samtools', 'view', '-b', bam_file, region, '-o', output_bam, '-@', str(20)]
    
    # 执行命令
    subprocess.run(cmd)
    print(f"BAM file created: {output_bam}")

# 使用samtools merge合并两个BAM文件
def merge_bam(output_bam, bam_part1, bam_part2):
    cmd = ['samtools', 'merge', output_bam, bam_part1, bam_part2]
    subprocess.run(cmd)
    print(f"Merged BAM file created: {output_bam}")

# 创建VCF文件
def create_vcf(row, unique_id, vcf_save_dir):
    header = [
        '##fileformat=VCFv4.2',
        '##source=Custom',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + row['Sample_ID']
    ]

    chrom_start = row['Chr_Start'].split(':')[0]
    pos_start = int(row['Chr_Start'].split(':')[1])
    chrom_end = row['Chr_End'].split(':')[0]
    pos_end = int(row['Chr_End'].split(':')[1])

    id_1 = f'{unique_id}_1'
    id_2 = f'{unique_id}_2'
    ref = 'N'
    alt_1 = f'[{chrom_end}:{pos_end}[N'
    alt_2 = f'N[{chrom_start}:{pos_start}]'
    qual = '.'
    filter_ = 'PASS'
    info_1 = (f'SVTYPE={row["Mutation_Type"]};STRANDS=--:503;EVENT={unique_id};MATEID={id_2};'
              f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={row["Frequency(%)"]};'
              f'PE=0;SR={row["Mean_Depth_Dedup"]};END={pos_end}')
    info_2 = (f'SVTYPE={row["Mutation_Type"]};STRANDS=--:503;EVENT={unique_id};MATEID={id_1};'
              f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={row["Frequency(%)"]};'
              f'PE=0;SR={row["Mean_Depth_Dedup"]};END={pos_start};SECONDARY')
    format_ = 'GT:SU:PE:SR'
    sample_data = f'./.:{row["Frequency(%)"]}:0:{row["Mean_Depth_Dedup"]}'
    
    vcf_data = [
        f'{chrom_start}\t{pos_start}\t{id_1}\t{ref}\t{alt_1}\t{qual}\t{filter_}\t{info_1}\t{format_}\t{sample_data}',
        f'{chrom_end}\t{pos_end}\t{id_2}\t{ref}\t{alt_2}\t{qual}\t{filter_}\t{info_2}\t{format_}\t{sample_data}'
    ]
    
    vcf_content = '\n'.join(header + vcf_data)
    with open(os.path.join(vcf_save_dir, f'{unique_id}.vcf'), 'w') as f:
        f.write(vcf_content)

# 处理VCF文件和裁剪BAM文件
def process_vcf(file_path, bam_dir, vcf_save_dir, bam_save_dir, encoding='latin1'):
    df = read_csv(file_path, encoding)
    
    # 确保保存VCF和BAM文件的目录存在
    ensure_directory_exists(vcf_save_dir)
    ensure_directory_exists(bam_save_dir)
    
    # 为每个样本创建编号或使用突变基因作为尾缀
    df['Unique_ID'] = df.groupby('Sample_ID').cumcount() + 1
    df['Unique_Suffix'] = df['Mutation_Base'] + '_' + df['Unique_ID'].astype(str)
    df['Unique_ID'] = df['Sample_ID'] + '_' + df['Unique_Suffix']
    
    # 保存新的CSV文件
    numbered_csv_path = os.path.join(vcf_save_dir, 'numbered_variants.csv')
    df.to_csv(numbered_csv_path, index=False)
    
    # 为每个变异生成对应的VCF文件，并裁剪BAM文件
    for idx, row in df.iterrows():
        unique_id = row['Unique_ID']
        
        # 生成VCF文件，保存到指定的VCF目录
        create_vcf(row, unique_id, vcf_save_dir)
        
        # 裁剪BAM文件, 针对起始和结束位置裁剪，保存到指定的BAM目录
        chrom_start = row['Chr_Start'].split(':')[0]
        pos_start = int(row['Chr_Start'].split(':')[1])
        chrom_end = row['Chr_End'].split(':')[0]
        pos_end = int(row['Chr_End'].split(':')[1])
        
        # 起始位置前1000000和后10000
        output_bam_start = os.path.join(bam_save_dir, f"{unique_id}_start.bam")
        clip_bam(row['Sample_ID'], chrom_start, max(0, pos_start - 10000), pos_start + 1000000, bam_dir, output_bam_start)
        
        # 结束位置前1000000和后10000
        output_bam_end = os.path.join(bam_save_dir, f"{unique_id}_end.bam")
        clip_bam(row['Sample_ID'], chrom_end, max(0, pos_end - 1000000), pos_end + 10000, bam_dir, output_bam_end)
        
        # 合并两个BAM文件，保存到指定的BAM目录
        merged_bam = os.path.join(bam_save_dir, f"{unique_id}.bam")
        merge_bam(merged_bam, output_bam_start, output_bam_end)
        
        # 删除中间的部分文件
        os.remove(output_bam_start)
        os.remove(output_bam_end)

# 使用示例
file_path = "/data/home/std_12/ShiHe/sarcoma.csv"
bam_dir = "/data/home/std_12/ShiHe/sarcoma_BAM"   # BAM文件的目录
vcf_save_dir = "/data/home/std_12/ShiHe/sarcomavcf"   # VCF文件保存目录
bam_save_dir = "/data/home/std_12/ShiHe/sarcomabam"   # BAM文件保存目录

process_vcf(file_path, bam_dir, vcf_save_dir, bam_save_dir)
print('wancheng')