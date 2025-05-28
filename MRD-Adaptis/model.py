import torch
import torch.nn as nn
import csv
import argparse
import os
import pysam
from tqdm import tqdm
import concurrent.futures
import subprocess
import argparse
import time
import pandas as pd
def convert_bam_to_sam(input_bam, output_sam):
    try:
        # 构建samtools命令
        samtools_command = ["samtools", "view", "-h", "-o", output_sam, input_bam]

        # 执行samtools命令
        subprocess.run(samtools_command, check=True)

        print(f"Conversion from BAM to SAM successful. SAM file saved as {output_sam}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running samtools: {e}")


def is_number(s):
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
    sss = 0
    for iii in cigar_arr:
        if iii[-1] == "M" or iii[-1] == "D":
            sss = sss + int(iii[:-1])
    return sss

def process_sam_file(file_name):
    this_samFile_arr1 = []
    this_samFile_readID_maxMapQul_dic = {}

    with open(file_name, 'r') as this_samFile_obj:
        for oneSamLine in this_samFile_obj.readlines():
            if not oneSamLine.startswith('@'):
                oneSamLine_array = oneSamLine.split('\t')
                if oneSamLine_array[5] != "*" and int(oneSamLine_array[4]) >= 30 and "H" not in oneSamLine_array[5]:
                    if oneSamLine_array[0] in this_samFile_readID_maxMapQul_dic:
                        if int(oneSamLine_array[4]) > this_samFile_readID_maxMapQul_dic[oneSamLine_array[0]]:
                            this_samFile_readID_maxMapQul_dic[oneSamLine_array[0]] = int(oneSamLine_array[4])
                    else:
                        this_samFile_readID_maxMapQul_dic[oneSamLine_array[0]] = int(oneSamLine_array[4])
                    this_samFile_arr1.append(oneSamLine_array)
    
   # print("第一步完成")

    this_samFile_arr2 = []
    for oneSamLine1 in this_samFile_arr1:
        if this_samFile_readID_maxMapQul_dic[oneSamLine1[0]] == int(oneSamLine1[4]):
            this_samFile_arr2.append(oneSamLine1)
    
   # print("第二步完成")
    
    return this_samFile_arr2

def process_structural_variations(this_samFile_arr2, chromosomesCodingDic):
    sv_begin_arr = []
    sv_end_arr = []
    read_len_arr = []

    for oneSamLine2 in this_samFile_arr2:
        chromosomes = oneSamLine2[2]
        if "_" in chromosomes:
            chromosomes = chromosomes.split("_")[0]
        ini_mapPos = int(oneSamLine2[3])
        mapPos = chromosomesCodingDic[chromosomes] + ini_mapPos 
        mapCigar = oneSamLine2[5]
        cigar_arr = []
        cigar_arr_one = ""
        for one in mapCigar:
            cigar_arr_one = cigar_arr_one + one
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
    threshold_value = 110000
    if end_pos > begin_pos and (end_pos - begin_pos) < threshold_value:
        return True
    else:
        return False
def find_similar_svs(sv_begin_arr, sv_end_arr):
    sv_result = []
    for one_sv_begin in sv_begin_arr:
        for one_sv_end in sv_end_arr:
            if same_sv(one_sv_begin, one_sv_end):
                sv_result.append([one_sv_begin, one_sv_end])
    return sv_result
def read_repeat_file(iniRepeatFileDir, chromosomesCodingDic):
    try:
        with open(iniRepeatFileDir, 'r') as iniRepeatObject:
            iniRepeatFile = iniRepeatObject.read()
    except Exception as e:
    #    print(f"An error occurred while reading the repeat file: {str(e)}")
        iniRepeatFile = ""

    iniRepeatFile_rows = iniRepeatFile.split('\n')
    iniRepeatFile_rows.pop()
    repeat_arr = []
    for one_iniRepeatFile_row in iniRepeatFile_rows:
        one_iniRepeatFile_row_arr = one_iniRepeatFile_row.split('\t')
        try:
            repeat_arr.append([chromosomesCodingDic[one_iniRepeatFile_row_arr[5]] + int(one_iniRepeatFile_row_arr[6]), chromosomesCodingDic[one_iniRepeatFile_row_arr[5]] + int(one_iniRepeatFile_row_arr[7])])
        except:
            continue

    return repeat_arr

def calculate_repeat_percentage(sv_result, repeat_arr):
    svInRepeat_account = 0
    for one_sv in sv_result:
        for one_repeat in repeat_arr:
            if same_sv(one_sv[0], one_sv[1]):
                svInRepeat_account = svInRepeat_account + 1
                break
    
#    repeatPercent = float(svInRepeat_account) / float(len(sv_result))
    if len(sv_result) == 0:
        repeatPercent = 0.0  # 或者其他你认为合适的默认值
    else:
        repeatPercent = float(svInRepeat_account) / float(len(sv_result))

    return repeatPercent
def classify_sv_length(sv_result):
    shortSV_account = 0
    middleSV_account = 0
    longSV_account = 0
    for one_sv in sv_result:
        if (one_sv[1] - one_sv[0]) <= 200:
            shortSV_account = shortSV_account + 1
        elif (one_sv[1] - one_sv[0]) > 1000:
            longSV_account = longSV_account + 1
        else:
            middleSV_account = middleSV_account + 1

    total_sv = len(sv_result)
#    total_sv = len(sv_result)

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
    if len(read_len_arr) == 0:
        return 0
    else:
        return float(sum(read_len_arr)) / float(len(read_len_arr))

def calculate_small_gap_metrics(sv_result,readLen):
    gap_arr = []
    before_sv_begin = None

    for one_sv in sv_result:
        if before_sv_begin is None:
            before_sv_begin = one_sv[0]
        else:
            now_sv_begin = one_sv[0]
            if str(now_sv_begin)[0:-9] == str(before_sv_begin)[0:-9]:
                gap_arr.append(now_sv_begin - before_sv_begin)
            before_sv_begin = now_sv_begin

    smallGap_account = 0
    smallGap_arr = []

    for one_gap in gap_arr:
        if one_gap <= min(readLen, 20000):
            smallGap_account = smallGap_account + 1
            smallGap_arr.append(one_gap)

    if len(smallGap_arr) == 0:
        RMB = 1.0
    else:
        smallGap_ave = float(sum(smallGap_arr)) / float(len(smallGap_arr))
        if smallGap_ave == 0.0:
            RMB = 1.0  # 或者您认为合适的非零值
        else:
            RMB = float(readLen) / float(smallGap_ave)

#    HMDP = float(len(smallGap_arr)) / float(len(gap_arr))
    if len(gap_arr) == 0:
        HMDP =1.0  # 或者其他你认为合适的默认值
    else:
        HMDP = float(len(smallGap_arr)) / float(len(gap_arr))


    return smallGap_account, RMB, HMDP

def calculate_average_depth(depth_file):
    total_depth = 0
    num_positions = 0
    with open(depth_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                depth = int(parts[2])
            except ValueError:
                continue
            total_depth += depth
            num_positions += 1
    if num_positions == 0:
        return 0
    return total_depth / num_positions

def calculate_average_depth_from_bam(bam_file):
    base_name = os.path.basename(bam_file)  # 提取文件名，例如 "1.bam"
    tmp_file = base_name + "_tmp_depth.txt"  # 生成临时文件名，例如 "1.bam_tmp_depth.txt"

    # 生成深度信息
    subprocess.run(["samtools", "depth", bam_file], stdout=open(tmp_file, 'w'))

    # 从深度文件计算平均深度
    average_depth = calculate_average_depth(tmp_file)

    # 删除临时文件（如果需要的话）
    os.remove(tmp_file)

    return average_depth


def process(input_bam_file):
    # 将BAM文件转换为SAM文件
    chromosomesCodingDic = {"chr1": 1000000000, "chr2": 2000000000, "chr3": 3000000000, "chr4": 4000000000, "chr5": 5000000000, "chr6": 6000000000, "chr7": 7000000000, "chr8": 8000000000, "chr9": 9000000000, "chr10": 10000000000, "chr11": 11000000000, "chr12": 12000000000, "chr13": 13000000000, "chr14": 14000000000, "chr15": 15000000000, "chr16": 16000000000, "chr17": 17000000000, "chr18": 18000000000, "chr19": 19000000000, "chr20": 20000000000, "chr21": 21000000000, "chr22": 22000000000, "chrX": 23000000000, "chrY": 24000000000, "chrM": 25000000000, "1": 1000000000, "2": 2000000000, "3": 3000000000, "4": 4000000000, "5": 5000000000, "6": 6000000000, "7": 7000000000, "8": 8000000000, "9": 9000000000, "10": 10000000000, "11": 11000000000, "12": 12000000000, "13": 13000000000, "14": 14000000000, "15": 15000000000, "16": 16000000000, "17": 17000000000, "18": 18000000000, "19": 19000000000, "20": 20000000000, "21": 21000000000, "22": 22000000000, "X": 23000000000, "Y": 24000000000, "M": 25000000000}
    iniRepeatFileDir = "rmsk.txt"
    output_sam_file = input_bam_file.replace(".bam", ".sam")
    convert_bam_to_sam(input_bam_file, output_sam_file)

    # 处理SAM文件
    this_samFile = output_sam_file
    result_arr2 = process_sam_file(this_samFile)
    sv_begin_arr, sv_end_arr, read_len_arr = process_structural_variations(result_arr2, chromosomesCodingDic)
    sv_result = find_similar_svs(sv_begin_arr, sv_end_arr)
    repeat_arr = read_repeat_file(iniRepeatFileDir, chromosomesCodingDic)
    repeatPercent = calculate_repeat_percentage(sv_result, repeat_arr)
    shortSV_percent, middleSV_percent, longSV_percent = classify_sv_length(sv_result)
    readLen = calculate_read_length(read_len_arr)
    smallGap_account, RMB, HMDP = calculate_small_gap_metrics(sv_result,readLen)
    average_depth = calculate_average_depth_from_bam(input_bam_file)

    # 定义要保存的数据
    data = {
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
    # 从字典中提取值并转换为列表
    data_values = list(data.values())

    # 将列表转换为PyTorch张量
    data_tensor = torch.tensor(data_values, dtype=torch.float32)

    return data_tensor

start_time = time.time()

n_features = 9  # 示例特征数量
n_targets = 8    # 示例目标数量
# 确保有模型类的定义
class MultiTargetRegression(nn.Module):
    def __init__(self):
        super(MultiTargetRegression, self).__init__()
        self.fc1 = nn.Linear(n_features, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, n_targets)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x

# 加载模型

model = MultiTargetRegression()
model.load_state_dict(torch.load('multi_target_regression_model1.pth'))
model.eval()  # 设置为评估模式

# 使用模型进行预测
# 示例输入（你需要用实际的数据替换这部分）
# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description="Process ")
# 添加 'input_file' 命令行参数
parser.add_argument('input_file', type=str, help="Path to the input file")
# 解析命令行参数
args = parser.parse_args()
example_input = process(args.input_file)  # 有一个样本，n_features特征
# print(example_input)
# print(example_input)
prediction = model(example_input).tolist()
# 创建一个字典，每个维度对应一个名称
# 创建包含列表的字典
prediction_dict = {
    "msw": [prediction[0]],
    "tt": [prediction[1]],
    "back_distance": [prediction[2]],
    "min_mapping_threshold": [prediction[3]],
    "min_clip": [prediction[4]],
    "read_length": [prediction[5]],
    "min_non_overlap": [prediction[6]],
    "discordant_z": [prediction[7]]
}
prediction_df = pd.DataFrame(prediction_dict)
print("The recommended parameter strategy is:")
print(prediction_df)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"运行时间: {elapsed_time} 秒")