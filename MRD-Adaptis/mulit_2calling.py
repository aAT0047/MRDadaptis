#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import subprocess
import re
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import csv
import sys
import time
import random
from evaf1 import *
from config import *
import re
import shutil
# Define the global variables
# bam_dir = "/home/cloudam/simulat_2/sample_bam"
# vcf_dir = "/home/cloudam/simulat_2/sample_vcf"
# 从CSV文件读取前305个样本的数据

# 定义参数范围
param_bounds = {
    "w": (1000, 10000),
    "msw": (1, 6),
    "tt": (0, 50),
    "back_distance": (0, 150),
    "min_mapping_threshold": (0, 100),
    "min_clip": (10, 500),
    "read_length": (50, 200),
    "min_non_overlap": (50, 200),
    "discordant_z": (2,9)
}



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
        id_match = re.search(r'EVENT=([^;]+)', info)
        su_match = re.search(r'SU=(\d+)', info)
        if id_match and su_match:
            event_id = id_match.group(1)
            su = int(su_match.group(1))
            if event_id not in variants:
                variants[event_id] = []
            variants[event_id].append((line, su))

    # 按支持证据数量排序，选择支持证据最多的变异对
    best_variant_pair = []
    highest_su = 0
    for event_id, var_list in variants.items():
        total_su = sum([su for _, su in var_list])
        if total_su > highest_su:
            highest_su = total_su
            best_variant_pair = var_list

    # 生成过滤后的临时 VCF 文件
    with open(temp_vcf_file, 'w') as f:
        for line in header:
            f.write(line)
        for variant, _ in best_variant_pair:
            f.write(variant)

    # 用临时 VCF 文件覆盖原始 VCF 文件
    shutil.move(temp_vcf_file, vcf_file)


def read_and_clean_csv(filename):
    print('开始清理分布小于10的数据')

    # 读取csv文件到DataFrame
    df = pd.read_csv(filename)
    # 处理莫名其妙失踪的一个bam
    # df = df[df['Sample'] == 'TE23B60067-D2ABAP2GXF1-L000_ROS1_1']
    # 删除含有空值的行
    df_cleaned = df.dropna()
    print(df_cleaned.columns)
    
    # 将mean和stdev列转换为浮点数
    # df_cleaned['Mean'] = df_cleaned['Mean'].astype(float)
    # df_cleaned['Stdev'] = df_cleaned['Stdev'].astype(float)
    df_cleaned.loc[:, 'Mean'] = df_cleaned['Mean'].astype(float)
    df_cleaned.loc[:, 'Stdev'] = df_cleaned['Stdev'].astype(float)

    # 过滤掉 'Sample' 列值为 100 的行 (uncomment this line if needed)
    # df_cleaned = df_cleaned[df_cleaned['Sample'] == 100]

    # 保存处理后的数据回csv文件
    df_cleaned.to_csv(filename, index=False)
    print('清理数据结束')

    # 将每行转换为元组并存入列表
    samples = list(df_cleaned.itertuples(index=False, name=None))
    # print(samples)
    return samples

def execute_command(cmd):
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        return result.decode('utf-8').strip()
    except subprocess.CalledProcessError as e:
        return e.output.decode('utf-8').strip()

def init_tqdm(pbar):
    def tqdm_callback(*a, **k):
        pbar.update()
    return tqdm_callback



def process_bam_files(sample_name):
    # bam_dir  = "/home/cloudam/simulat_2/sample_bam"
    
    # 提取discordant配对末端比对数据
    bam_file = os.path.join(bam_dir, f"{sample_name}.bam")
    
    # 索引初始的BAM文件
    subprocess.call(["samtools", "index", bam_file])
    
    discordants_unsorted_file = os.path.join(bam_dir, f"{sample_name}.discordants.unsorted.bam")
    cmd = f"samtools view -b -F 1294 {bam_file} > {discordants_unsorted_file}"
    subprocess.call(cmd, shell=True)

    # Extract the split-read alignments
    splitters_unsorted_file = os.path.join(bam_dir, f"{sample_name}.splitters.unsorted.bam")
    cmd = (
        f"samtools view -h {bam_file} | "
        f"/data/home/std_12/lumpy-sv-0.3.1/scripts/extractSplitReads_BwaMem -i stdin | "
        f"samtools view -Sb - > {splitters_unsorted_file}"
    )
    subprocess.call(cmd, shell=True)

    # 排序discordant和split-read数据
    discordants_file = os.path.join(bam_dir, f"{sample_name}.discordants.bam")
    subprocess.call(["samtools", "sort", "-o", discordants_file, discordants_unsorted_file])
    
    splitters_file = os.path.join(bam_dir, f"{sample_name}.splitters.bam")
    subprocess.call(["samtools", "sort", "-o", splitters_file, splitters_unsorted_file])
    
    # 索引sorted BAM文件
    subprocess.call(["samtools", "index", discordants_file])
    subprocess.call(["samtools", "index", splitters_file])

    # 清理未排序的文件
    os.remove(discordants_unsorted_file)
    os.remove(splitters_unsorted_file)


def generate_random_params(param_bounds):
    """生成随机参数列表"""
    return {
        "w": random.randint(*param_bounds["w"]),
        "msw": random.randint(*param_bounds["msw"]),
        "tt": random.randint(*param_bounds["tt"]),
        "back_distance": random.randint(*param_bounds["back_distance"]),
        "min_mapping_threshold": random.randint(*param_bounds["min_mapping_threshold"]),
        "min_clip": random.randint(*param_bounds["min_clip"]),
        "read_length": random.randint(*param_bounds["read_length"]),
        "min_non_overlap": random.randint(*param_bounds["min_non_overlap"]),
        "discordant_z": random.randint(*param_bounds["discordant_z"]),
    }

def write_log(log_file, message):
    """写入日志文件"""
    with open(log_file, 'a') as f:
        f.write(message + '\n')

def move_file(src, dest):
    """移动文件"""
    mv_cmd = f'mv {src} {dest}'
    subprocess.call(mv_cmd, shell=True)



def process_sample(data):
    sample_name, mean, stdev = data
    sample_name = str(sample_name)
    output_vcf_path = os.path.join(vcf_dir, f"{sample_name}.vcf")

    log_file = os.path.join(log_folder, f"{sample_name}_log.txt")

    write_log(log_file, f"开始处理样本 {sample_name}")

    if not os.path.exists(csv_path):
        os.makedirs(csv_path, exist_ok=True)

    initialization_csv = os.path.join(csv_path, f"{sample_name}_params.csv")

    with open(initialization_csv, 'w') as f:
        f.write("sample_name,w,msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,mean, stdev,TP,FP,FN,f1_score,precision,recall\n")

    random_params_list = [generate_random_params(param_bounds) for _ in range(1999)]
    # 指定的参数集
    specific_params = {
        "w": 1000000,
        "msw": 4,
        "tt": 0,
        "back_distance": 10,
        "min_mapping_threshold": 20,
        "min_clip": 50,
        "read_length": 150,
        "min_non_overlap": 150,
        "discordant_z": 5
    }

    # 将指定的参数集添加到列表中
    random_params_list.append(specific_params)

    current_iteration = 0
    total_iterations = len(random_params_list)

    new_folder = os.path.join(bam_dir, sample_name)
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)

    os.chdir(new_folder)

    for params in random_params_list:
        w = params["w"]
        msw = params["msw"]
        tt = params["tt"]
        back_distance = params["back_distance"]
        min_mapping_threshold = params["min_mapping_threshold"]
        min_clip = params["min_clip"]
        read_length = params["read_length"]
        min_non_overlap = params["min_non_overlap"]
        discordant_z = params["discordant_z"]

        output_vcf = os.path.join(bam_dir, f"{sample_name}_variants_w{w}_msw{msw}_tt{tt}_bd{back_distance}_mmt{min_mapping_threshold}_mc{min_clip}_rl{read_length}_mno{min_non_overlap}_dz{discordant_z}.vcf")

        log_message = f"正在运行参数设置：w={w}, msw={msw}, tt={tt}, back_distance={back_distance}, min_mapping_threshold={min_mapping_threshold}, min_clip={min_clip}, mean={mean}, stdev={stdev}, read_length={read_length}, min_non_overlap={min_non_overlap}, discordant_z={discordant_z}"
        write_log(log_file, log_message)

        cmd = [
            "lumpy",
            "-w", str(w),
            "-msw", str(msw),
            "-tt", str(tt),
            "-pe", f"id:{sample_name},read_group:readgroup1,bam_file:{bam_dir}/{sample_name}.discordants.bam,histo_file:{bam_dir}/{sample_name}.lib1.histo,mean:{mean},stdev:{stdev},read_length:{read_length},min_non_overlap:{min_non_overlap},discordant_z:{discordant_z},back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold}",
            "-sr", f"id:{sample_name},bam_file:{bam_dir}/{sample_name}.splitters.bam,back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold}",
            ">", output_vcf
        ]

        # 打开日志文件
        with open(log_file, 'w') as log:
            # 运行命令并将输出重定向到日志文件
            process = subprocess.Popen(" ".join(cmd), shell=True, stdout=log, stderr=log)
            process.communicate()
        filter_vcf_by_evidence(output_vcf)
        TP, FP, FN, f1_score, precision, recall = evamain(output_vcf_path, output_vcf)

        move_file(output_vcf, last_vcf_dir)

        with open(initialization_csv, 'a') as f:
            f.write(f"{sample_name},{w},{msw},{tt},{back_distance},{min_mapping_threshold},{min_clip},{read_length},{min_non_overlap},{discordant_z},{mean},{stdev},{TP},{FP},{FN},{f1_score},{precision},{recall}\n")

        current_iteration += 1
        sys.stdout.write(f"\r进度: {current_iteration}/{total_iterations} ({(float(current_iteration) / total_iterations) * 100:.2f}%)")
        sys.stdout.flush()

    print(f"Sample {sample_name} processed.")

def mulit_2lumpymain():

    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    sample_names = read_and_clean_csv(csv_file)
    sample_datas = [row[0] for row in sample_names]  # 提取每一行的第一个元素，即样本名

    print('准备discordant和split-read数据')
    with tqdm(total=len(sample_datas), desc="Processing BAM files") as pbar:
        with Pool(processes=cpu_count()) as pool:
            for _ in pool.imap_unordered(process_bam_files, sample_datas):
                pbar.update()  # 更新进度条
    print('discordant和split-read数据完成')

    print('开始运行LUMPY')
    with tqdm(total=len(sample_names), desc="Processing samples") as pbar:
        with Pool(processes=cpu_count()) as pool:
            for _ in pool.imap(process_sample, sample_names):
                pbar.update()  # 更新进度条
    print("All samples processed.")