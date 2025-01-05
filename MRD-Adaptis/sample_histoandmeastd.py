import os
import subprocess
import tempfile
import config
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd
import re

def load_samtools():
    # 加载 SAMtools 模块
    subprocess.check_call("module load SAMtools/1.15.1-GCC-11.2.0", shell=True)

def create_log_file(log_file, summary_csv_file):
    # 创建日志文件
    with open(log_file, 'w') as f:
        f.write("Processing log\n")
        f.write(f"Logs will be written to {log_file}\n")
        f.write(f"Summary CSV file will be written to {summary_csv_file}\n")

def get_bam_files(input_dir):
    # 获取所有BAM文件，并限制只获取前100个文件
    return [f for f in os.listdir(input_dir) if f.endswith('.bam')]

def process_sample(bam_file, input_dir, output_dir, log_file):
    sample_name = os.path.splitext(bam_file)[0]
    csv_file = os.path.join(output_dir, f"{sample_name}.csv")
    histo_file = os.path.join(input_dir, f"{sample_name}.lib1.histo")
    # 跳过包含指定后缀的文件
    if any(ext in bam_file for ext in [".discordants.bam", ".splitters.bam", 
                                       ".discordants.unsorted", ".splitters.unsorted", 
                                       ".discordants.sorted", ".splitters.sorted"]):
        with open(log_file, 'a') as f:
            f.write(f"Skipping {bam_file}\n")
        return

    # 使用 tempfile 模块生成临时文件
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        try:
            temp_file_path = temp_file.name
            # 生成经验插入大小统计数据
            subprocess.check_call(f"samtools view {os.path.join(input_dir, bam_file)} > {temp_file_path}", shell=True)
            
            output = subprocess.check_output(
                f"tail -n+10000 {temp_file_path} | /data/home/std_12/lumpy-sv-0.3.1/scripts/pairend_distro.py -r 151 -X 4 -N 1000 -o {histo_file} 2>/dev/null", 
                shell=True
            ).decode('utf-8')

            mean = "NA"
            stdev = "NA"
            
            for line in output.split('\n'):
                if 'mean:' in line:
                    mean = line.split('mean:')[1].split()[0]
                if 'stdev:' in line:
                    stdev = line.split('stdev:')[1].split()[0]

            stdev = re.sub(r'stdev:', '', stdev)

            if mean == "NA" or stdev == "NA":
                with open(log_file, 'a') as f:
                    f.write(f"Failed to get mean or stdev for {sample_name}\n")

            # 写入单独的CSV文件
            with open(csv_file, 'w') as f:
                f.write("Sample,Mean,Stdev\n")
                f.write(f"{sample_name},{mean},{stdev}\n")

            with open(log_file, 'a') as f:
                f.write(f"Sample {sample_name} processed.\n")
        
        except subprocess.CalledProcessError as e:
            with open(log_file, 'a') as f:
                f.write(f"Error processing {bam_file}: {str(e)}\n")
        finally:
            os.remove(temp_file_path)

def merge_csv_files(output_dir, summary_csv_file):
    # 合并所有单独的CSV文件为一个汇总文件
    all_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('.csv')]
    combined_csv = pd.concat([pd.read_csv(f) for f in all_files])
    combined_csv.to_csv(summary_csv_file, index=False)

    # 删除临时CSV文件
    for f in all_files:
        os.remove(f)

def sample_histoandmeastdmain():
    # 从配置文件中获取全局变量
    input_dir = config.input_dir
    output_dir = config.output_dir
    summary_csv_file = config.csv_file
    log_file = config.log_file

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # load_samtools()
    create_log_file(log_file, summary_csv_file)
    bam_files = get_bam_files(input_dir)

    # 使用多线程处理样本并显示进度条
    with ThreadPoolExecutor(max_workers=32) as executor:
        futures = {executor.submit(process_sample, bam_file, input_dir, output_dir, log_file): bam_file for bam_file in bam_files}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing BAM files"):
            future.result()

    merge_csv_files(output_dir, summary_csv_file)

    print(f"All samples processed. Summary CSV file is located at {summary_csv_file}")


