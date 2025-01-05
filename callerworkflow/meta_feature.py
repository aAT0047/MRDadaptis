import csv
import os
import pysam
from tqdm import tqdm
import concurrent.futures
import subprocess
import pandas as pd

# 定义染色体编码字典
chromosomesCodingDic = {
    "chr1": 1000000000, "chr2": 2000000000, "chr3": 3000000000, "chr4": 4000000000,
    "chr5": 5000000000, "chr6": 6000000000, "chr7": 7000000000, "chr8": 8000000000,
    "chr9": 9000000000, "chr10": 10000000000, "chr11": 11000000000, "chr12": 12000000000,
    "chr13": 13000000000, "chr14": 14000000000, "chr15": 15000000000, "chr16": 16000000000,
    "chr17": 17000000000, "chr18": 18000000000, "chr19": 19000000000, "chr20": 20000000000,
    "chr21": 21000000000, "chr22": 22000000000, "chrX": 23000000000, "chrY": 24000000000,
    "chrM": 25000000000, "1": 1000000000, "2": 2000000000, "3": 3000000000, "4": 4000000000,
    "5": 5000000000, "6": 6000000000, "7": 7000000000, "8": 8000000000, "9": 9000000000,
    "10": 10000000000, "11": 11000000000, "12": 12000000000, "13": 13000000000,
    "14": 14000000000, "15": 15000000000, "16": 16000000000, "17": 17000000000,
    "18": 18000000000, "19": 19000000000, "20": 20000000000, "21": 21000000000,
    "22": 22000000000, "X": 23000000000, "Y": 24000000000, "M": 25000000000
}

# 定义配置集列表
configurations = [
    # {
    #     "name": "IGbam",
    #     "base_dir": "/data/home/std_12/ShiHe",
    #     "input_dir": "/data/home/std_12/ShiHe/IGbam",
    #     "output_base_dir": "/data/home/std_12/ShiHe/IGbam_csv"
    # },
    # {
    #     "name": "ALK_RET_ROS1bam",
    #     "base_dir": "/data/home/std_12/ShiHe",
    #     "input_dir": "/data/home/std_12/ShiHe/ALK-RET-ROS1bam",
    #     "output_base_dir": "/data/home/std_12/ShiHe/ALK_RET_ROS1bam_csv"
    # },
    {
        "name": "sarcomabam",
        "base_dir": "/data/home/std_12/ShiHe",
        "input_dir": "/data/home/std_12/ShiHe/sarcomabam",
        "output_base_dir": "/data/home/std_12/ShiHe/sarcomabam_csv"
    }
]

def convert_bam_to_sam(input_bam, output_sam):
    """
    将BAM文件转换为SAM文件
    """
    try:
        # 构建samtools命令
        samtools_command = ["samtools", "view", "-h", "-o", output_sam, input_bam]

        # 执行samtools命令
        subprocess.run(samtools_command, check=True)

        print(f"成功将BAM转换为SAM文件，保存为 {output_sam}")
    except subprocess.CalledProcessError as e:
        print(f"运行samtools时发生错误: {e}")

def is_number(s):
    """
    判断字符串是否为数字
    """
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def ref_span(cigar_arr):
    """
    计算CIGAR字符串对应的参考序列跨度
    """
    sss = 0
    for iii in cigar_arr:
        if iii[-1] == "M" or iii[-1] == "D":
            sss += int(iii[:-1])
    return sss

def process_sam_file(file_name):
    """
    处理SAM文件，过滤并提取所需信息
    """
    this_samFile_arr1 = []
    this_samFile_readID_maxMapQul_dic = {}

    with open(file_name, 'r') as this_samFile_obj:
        for oneSamLine in this_samFile_obj:
            if not oneSamLine.startswith('@'):
                oneSamLine_array = oneSamLine.strip().split('\t')
                if oneSamLine_array[5] != "*" and int(oneSamLine_array[4]) >= 30 and "H" not in oneSamLine_array[5]:
                    read_id = oneSamLine_array[0]
                    map_quality = int(oneSamLine_array[4])
                    if read_id in this_samFile_readID_maxMapQul_dic:
                        if map_quality > this_samFile_readID_maxMapQul_dic[read_id]:
                            this_samFile_readID_maxMapQul_dic[read_id] = map_quality
                    else:
                        this_samFile_readID_maxMapQul_dic[read_id] = map_quality
                    this_samFile_arr1.append(oneSamLine_array)

    # 第二步过滤
    this_samFile_arr2 = []
    for oneSamLine1 in this_samFile_arr1:
        if this_samFile_readID_maxMapQul_dic[oneSamLine1[0]] == int(oneSamLine1[4]):
            this_samFile_arr2.append(oneSamLine1)

    return this_samFile_arr2

def process_structural_variations(this_samFile_arr2):
    """
    处理结构变异，提取SV起始和结束位置以及读取长度
    """
    sv_begin_arr = []
    sv_end_arr = []
    read_len_arr = []

    for oneSamLine2 in this_samFile_arr2:
        chromosomes = oneSamLine2[2]
        if "_" in chromosomes:
            chromosomes = chromosomes.split("_")[0]
        ini_mapPos = int(oneSamLine2[3])
        mapPos = chromosomesCodingDic.get(chromosomes, 0) + ini_mapPos
        mapCigar = oneSamLine2[5]
        cigar_arr = []
        cigar_arr_one = ""
        for one in mapCigar:
            cigar_arr_one += one
            if not one.isdigit():
                cigar_arr.append(cigar_arr_one)
                cigar_arr_one = ""
        cigar_arr_first = cigar_arr[0]
        cigar_arr_last = cigar_arr[-1]
        if cigar_arr_first[-1] == "S" and int(cigar_arr_first[:-1]) >= 50:
            this_end_value = mapPos - 1
            end_flag = 0
            for one_sv_end in sv_end_arr:
                if abs(this_end_value - one_sv_end) < 51:
                    end_flag = 1
                    break
            if end_flag == 0:
                sv_end_arr.append(this_end_value)
                read_len = len(oneSamLine2[9])
                if read_len >= 100:
                    read_len_arr.append(read_len)
        if cigar_arr_last[-1] == "S" and int(cigar_arr_last[:-1]) >= 50:
            this_begin_value = mapPos + ref_span(cigar_arr)
            begin_flag = 0
            for one_sv_begin in sv_begin_arr:
                if abs(this_begin_value - one_sv_begin) < 51:
                    begin_flag = 1
                    break
            if begin_flag == 0:
                sv_begin_arr.append(this_begin_value)
                read_len = len(oneSamLine2[9])
                if read_len >= 100:
                    read_len_arr.append(read_len)

    sv_begin_arr = sorted(sv_begin_arr)
    sv_end_arr = sorted(sv_end_arr)

    return sv_begin_arr, sv_end_arr, read_len_arr

def same_sv(begin_pos, end_pos):
    """
    判断两个SV是否相同，阈值为110000
    """
    threshold_value = 110000
    if end_pos > begin_pos and (end_pos - begin_pos) < threshold_value:
        return True
    else:
        return False

def find_similar_svs(sv_begin_arr, sv_end_arr):
    """
    找到相似的SV对
    """
    sv_result = []
    for one_sv_begin in sv_begin_arr:
        for one_sv_end in sv_end_arr:
            if same_sv(one_sv_begin, one_sv_end):
                sv_result.append([one_sv_begin, one_sv_end])
    return sv_result

def read_repeat_file(iniRepeatFileDir):
    """
    读取重复序列文件，提取重复区域
    """
    try:
        with open(iniRepeatFileDir, 'r') as iniRepeatObject:
            iniRepeatFile = iniRepeatObject.read()
    except Exception as e:
        print(f"读取重复序列文件时发生错误: {str(e)}")
        iniRepeatFile = ""

    iniRepeatFile_rows = iniRepeatFile.strip().split('\n')
    repeat_arr = []
    for one_iniRepeatFile_row in iniRepeatFile_rows:
        one_iniRepeatFile_row_arr = one_iniRepeatFile_row.split('\t')
        try:
            chrom_code = chromosomesCodingDic.get(one_iniRepeatFile_row_arr[5], 0)
            repeat_arr.append([
                chrom_code + int(one_iniRepeatFile_row_arr[6]),
                chrom_code + int(one_iniRepeatFile_row_arr[7])
            ])
        except:
            continue

    return repeat_arr

def calculate_repeat_percentage(sv_result, repeat_arr):
    """
    计算SV在重复区域的比例
    """
    svInRepeat_account = 0
    for one_sv in sv_result:
        for one_repeat in repeat_arr:
            if same_sv(one_sv[0], one_sv[1]):
                svInRepeat_account += 1
                break

    if len(sv_result) == 0:
        repeatPercent = 0.0
    else:
        repeatPercent = float(svInRepeat_account) / float(len(sv_result))

    return repeatPercent

def classify_sv_length(sv_result):
    """
    分类SV长度，计算不同长度SV的比例
    """
    shortSV_account = 0
    middleSV_account = 0
    longSV_account = 0
    for one_sv in sv_result:
        sv_length = one_sv[1] - one_sv[0]
        if sv_length <= 200:
            shortSV_account += 1
        elif sv_length > 1000:
            longSV_account += 1
        else:
            middleSV_account += 1

    total_sv = len(sv_result)

    if total_sv == 0:
        shortSV_percent = 0.0
        middleSV_percent = 0.0
        longSV_percent = 0.0
    else:
        shortSV_percent = float(shortSV_account) / total_sv
        middleSV_percent = float(middleSV_account) / total_sv
        longSV_percent = float(longSV_account) / total_sv

    return shortSV_percent, middleSV_percent, longSV_percent

def calculate_read_length(read_len_arr):
    """
    计算平均读取长度
    """
    if len(read_len_arr) == 0:
        return 0
    else:
        return float(sum(read_len_arr)) / float(len(read_len_arr))

def calculate_small_gap_metrics(sv_result, readLen):
    """
    计算小间隔相关的指标
    """
    gap_arr = []
    before_sv_begin = None

    for one_sv in sv_result:
        if before_sv_begin is None:
            before_sv_begin = one_sv[0]
        else:
            now_sv_begin = one_sv[0]
            if str(now_sv_begin)[:-9] == str(before_sv_begin)[:-9]:
                gap_arr.append(now_sv_begin - before_sv_begin)
            before_sv_begin = now_sv_begin

    smallGap_account = 0
    smallGap_arr = []

    for one_gap in gap_arr:
        if one_gap <= min(readLen, 20000):
            smallGap_account += 1
            smallGap_arr.append(one_gap)

    if len(smallGap_arr) == 0:
        RMB = 1.0
    else:
        smallGap_ave = float(sum(smallGap_arr)) / float(len(smallGap_arr))
        if smallGap_ave == 0.0:
            RMB = 1.0
        else:
            RMB = float(readLen) / float(smallGap_ave)

    if len(gap_arr) == 0:
        HMDP = 1.0
    else:
        HMDP = float(len(smallGap_arr)) / float(len(gap_arr))

    return smallGap_account, RMB, HMDP

def calculate_average_depth_from_bam(bam_file):
    """
    计算BAM文件的平均测序深度
    """
    base_name = os.path.basename(bam_file)
    tmp_file = base_name + "_tmp_depth.txt"

    # 生成深度信息
    with open(tmp_file, 'w') as f:
        subprocess.run(["samtools", "depth", bam_file], stdout=f)

    # 计算平均深度
    total_depth = 0
    num_positions = 0
    with open(tmp_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                depth = int(parts[2])
            except ValueError:
                continue
            total_depth += depth
            num_positions += 1

    if num_positions == 0:
        average_depth = 0
    else:
        average_depth = total_depth / num_positions

    # 删除临时文件
    os.remove(tmp_file)

    return average_depth

def process_and_save_to_csv(input_bam_file, output_csv_file, sample_id):
    """
    处理单个BAM文件并将结果保存为CSV文件
    """
    # 将BAM文件转换为SAM文件
    output_sam_file = input_bam_file.replace(".bam", ".sam")
    convert_bam_to_sam(input_bam_file, output_sam_file)

    # 处理SAM文件
    result_arr2 = process_sam_file(output_sam_file)
    sv_begin_arr, sv_end_arr, read_len_arr = process_structural_variations(result_arr2)
    sv_result = find_similar_svs(sv_begin_arr, sv_end_arr)
    repeat_file_path = "/data/home/std_12/ShiHe/rmsk.txt"  # 修改为实际的rmsk.txt文件路径
    repeat_arr = read_repeat_file(repeat_file_path)
    repeatPercent = calculate_repeat_percentage(sv_result, repeat_arr)
    shortSV_percent, middleSV_percent, longSV_percent = classify_sv_length(sv_result)
    readLen = calculate_read_length(read_len_arr)
    smallGap_account, RMB, HMDP = calculate_small_gap_metrics(sv_result, readLen)
    average_depth = calculate_average_depth_from_bam(input_bam_file)

    # 定义要保存的数据
    data = {
        "Sample ID": sample_id,
        "Repeat Percentage": repeatPercent,
        "Short SV Percentage": shortSV_percent,
        "Middle SV Percentage": middleSV_percent,
        "Long SV Percentage": longSV_percent,
        "Average Read Length": readLen,
        "Small Gap Account": smallGap_account,
        "RMB": RMB,
        "HMDP": HMDP,
        "Average Depth": average_depth
    }

    # 写入CSV文件前，确保目录存在
    output_csv_dir = os.path.dirname(output_csv_file)
    os.makedirs(output_csv_dir, exist_ok=True)

    # 写入CSV文件
    with open(output_csv_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=data.keys())

        # 写入CSV文件的标题行
        writer.writeheader()

        # 写入数据行
        writer.writerow(data)

    print(f"样本 {sample_id} 的数据已保存到 {output_csv_file}")

def get_bam_files(input_dir):
    """
    获取指定目录下的所有BAM文件
    """
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    return bam_files

def main():
    for config in configurations:
        # 获取当前配置的名称，用于标识
        config_name = config["name"]
        print(f"\n正在处理配置集：{config_name}")

        # 定义要创建的输出文件夹路径
        output_data_dir = config["output_base_dir"]

        # 使用os.makedirs()创建输出文件夹，如果文件夹已存在则不会重复创建
        os.makedirs(output_data_dir, exist_ok=True)

        # 获取所有BAM文件
        bam_files = get_bam_files(config["input_dir"])
        # print(bam_files)
    
        # 创建用于存储单个CSV文件的目录
        csv_output_dir = os.path.join(output_data_dir, "csv_files")
        os.makedirs(csv_output_dir, exist_ok=True)

        # 使用多核并行处理多个BAM文件
        with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
            bam_paths = [os.path.join(config["input_dir"], bam_file) for bam_file in bam_files]
            csv_files = [os.path.join(csv_output_dir, bam_file.replace('.bam', '_data.csv')) for bam_file in bam_files]
            sample_ids = [bam_file.replace('.bam', '') for bam_file in bam_files]

            futures = []
            for bam_file, csv_file, sample_id in zip(bam_paths, csv_files, sample_ids):
                futures.append(executor.submit(process_and_save_to_csv, bam_file, csv_file, sample_id))

            # 使用tqdm显示进度条
            for _ in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc=f"处理{config_name}的BAM文件"):
                pass

        # 合并所有CSV文件
        merged_csv_file = os.path.join(output_data_dir, f"{config_name}_metafeature.csv")
        csv_files = [os.path.join(csv_output_dir, f.replace('.bam', '_data.csv')) for f in bam_files]
        merged_df = pd.read_csv(csv_files[0])
        # 检查第一个文件是否存在
        if os.path.exists(csv_files[0]):
            merged_df = pd.read_csv(csv_files[0])
        else:
            print(f"文件 {csv_files[0]} 不存在，跳过合并。")
            continue

        # 从第二个文件开始，读取并合并数据
        for csv_file in tqdm(csv_files[1:], desc=f"合并{config_name}的CSV文件"):
            if os.path.exists(csv_file):
                df = pd.read_csv(csv_file)
                merged_df = pd.concat([merged_df, df], ignore_index=True)
            else:
                print(f"文件 {csv_file} 不存在，跳过。")

        # 保存合并后的数据
        merged_df.to_csv(merged_csv_file, index=False)

        print(f"{config_name}的所有CSV文件已合并为 {merged_csv_file}")

        # # 删除单个CSV文件，只保留合并后的CSV
        # for csv_file in csv_files:
        #     os.remove(csv_file)
        # # 删除存储单个CSV文件的目录
        # os.rmdir(csv_output_dir)

    print("所有配置集处理完成。")

if __name__ == "__main__":
    main()
