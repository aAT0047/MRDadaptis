# -*- coding: UTF-8 -*
import re
import os
import random
import math
import time
import pandas as pd
import numpy as np
import random
import csv
from tqdm import tqdm

df =  pd.read_csv(r"C:\Users\1\Desktop\code_ppg\sv_sh\output_file_chr1.tsv", delimiter='\t')
def snv(Result = '',result=''):
    global SV_sim_begin
    for index, value in df['Variant_Type'].iteritems():
        if value == 'SNP':
            onesnv = 'onesnv'
            chromosome = df.iloc[index,df.columns.get_loc('Chromosome')]
            Slen =random.randint(30000, 50001)
            if df.iloc[index,df.columns.get_loc('Start_Position')] >=60000000 and df.iloc[index,df.columns.get_loc('End_Position')]  <=6100000000000000:
                chromosome_start = df.iloc[index,df.columns.get_loc('Start_Position')] 
                chromosome_end = df.iloc[index,df.columns.get_loc('End_Position')] 
                mutated_from_allele = df.iloc[index,df.columns.get_loc('Tumor_Seq_Allele1')]
                randomseed = random.random()
                mutated_to_allele = df.iloc[index,df.columns.get_loc('Tumor_Seq_Allele2')]
                list = [onesnv,str(chromosome),chromosome_start
                    ,mutated_from_allele, mutated_to_allele,1]
                converted_list = [str(item) for item in list]
                result = ','.join(converted_list)
                Result = Result + ":"+result
    return  Result[1:]
    
def getRef(pos_rnd):
    chr_rnd = "chr1"
    substitutition = {
        'A': ['T', 'C', 'G'],
        'T': ['A', 'C', 'G'],
        'C': ['T', 'A', 'G'],
        'G': ['T', 'C', 'A'],
        'a': ['T', 'C', 'G'],
        't': ['A', 'C', 'G'],
        'c': ['T', 'A', 'G'],
        'g': ['T', 'C', 'A'],
    }
    hg19index_dict = {
        "chr1": [5, 249250621],
        "chr2": [254235645, 243199373],
        "chr3": [502299012, 198022430],
        "chr4": [704281897, 191154276],
        "chr5": [899259265, 180915260],
        "chr6": [1083792837, 171115067],
        "chr7": [1258330212, 159138663],
        "chr8": [1420651655, 146364022],
        "chr9": [1569942964, 141213431],
        "chr10": [1713980671, 135534747],
        "chr11": [1852226120, 135006516],
        "chr12": [1989932774, 133851895],
        "chr13": [2126461714, 115169878],
        "chr14": [2243934997, 107349540],
        "chr15": [2353431535, 102531392],
        "chr16": [2458013562, 90354753],
        "chr17": [2550175418, 81195210],
        "chr18": [2632994540, 78077248],
        "chr19": [2712633340, 59128983],
        "chr20": [2772944910, 63025520],
        "chr21": [2837230948, 48129895],
        "chr22": [2886323448, 51304566],
        "chrX": [2938654112, 155270560],
        "chrY": [3097030090, 59373566]
    }
    hg_fa_file = file("hg19.fa", 'r')
    if (pos_rnd % 50) == 0:
        hg_fa_file.seek(hg19index_dict[chr_rnd][0] + pos_rnd + (pos_rnd / 50) - 1, 0)
    else:
        hg_fa_file.seek(hg19index_dict[chr_rnd][0] + pos_rnd + (pos_rnd / 50), 0)
    ref = hg_fa_file.read(2)
    ref = re.sub('\s', '', ref)
    ref = ref[0:1]
    ref = ref.upper()
    if ref == "N":
        alt = "N"
    else:
        substitute_index = int(round(random.random() * 2))
        alt = substitutition[ref][substitute_index]
    return ref, alt

def get_random_arr():
    element1 = random.randrange(900000, 1000001)
    element2 = random.randrange(900000, 1000001)
    element3 = random.randrange(900000, 1000001)
    element4 = random.randrange(900000, 1000001)
    elementSUM = element1 + element2 + element3 + element4
    element1proportion = float(element1)/float(elementSUM)
    element12proportion = element1proportion + float(element2)/float(elementSUM)
    element123proportion = element12proportion + float(element3)/float(elementSUM)
    element1234proportion = element123proportion + float(element4)/float(elementSUM)
    random_arr = [element1proportion, element12proportion, element123proportion, element1234proportion]
    return random_arr
output_directory = r"C:\Users\1\Desktop\code_ppg\sv_sh\sh_files"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
for iturn in tqdm(range(1, 2), desc="Generating .sh files"):  # for numbers 1 to 100
    # choices = list(range(90000000, 140000001, 10000000))
    # SV_sim_begin = random.choice(choices)
    # Create a new sh file for each iteration
    shFileName = os.path.join(output_directory, f'A_stableCallerPaperSimFlowShell_{iturn}.sh')
    sampleID = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
    shFile = open(shFileName, 'wb')
    depth_random = random.randrange(1, 21)
    this_sample_depth = str(depth_random * 10)
    # this_sample_readLength = str(random.choice([50, 75, 100, 150, 200, 250, 300]))
    this_sample_readLength  = 150

    #模拟结构变异开始
    SV_sim_begin = 10000000
    tsml_arr = get_random_arr()       #tinny, short, middle, large长度SV的累加占比数组
    tsml_arr[1] = tsml_arr[1] + math.floor((0.9-tsml_arr[1])*10)/10
    tsml_arr[2] = tsml_arr[2] + math.floor((1-tsml_arr[2])*10)/10
    diti_arr = get_random_arr()       #deletion, insertion, tandemRepeat, inversion类型SV的累加占比数组
    SV_block = random.choice([1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000])   #随机生成放置一个sv的block大小，用来控制变异密度
    deletionShellCMD = ""
    insertionShellCMD = ""
    tandemRepeatShellCMD = ""
    inversionShellCMD = ""
    want_to_sim_SV_number = random.randrange(450, 550)
    for i in range(want_to_sim_SV_number):
        tsml_random = random.random()
        diti_random = random.random()
        this_sv_begin = random.randrange(SV_sim_begin + SV_block * i, SV_sim_begin + SV_block * (i + 1))
        if tsml_random <= tsml_arr[0]:
            this_sv_len = random.randrange(900000, 1000001)
        elif tsml_random <= tsml_arr[1]:
            this_sv_len = random.randrange(900000, 1000001)
        elif tsml_random <= tsml_arr[2]:
            this_sv_len = random.randrange(900000, 1000001)
        else:
            this_sv_len = random.randrange(900000, 1000001)
        this_sv_end = min(this_sv_begin + this_sv_len, SV_sim_begin + SV_block * (i + 1) - 1)
        this_sv_begin = str(this_sv_begin)
        this_sv_end = str(this_sv_end)
        
        if diti_random <= diti_arr[0]:
            sv_type = "deletion"
            ifDelins_random = random.random()
            if ifDelins_random < 0.97:
                deletionShellCMD = deletionShellCMD + "onedelete,chr1," + this_sv_begin + "," + this_sv_end + ",1,0:"
            else:
                deletionShellCMD = deletionShellCMD + "onedelete,chr1," + this_sv_begin + "," + this_sv_end + ",1,1:"
        elif diti_random <= diti_arr[1]:
            sv_type = "insertion"
            insertionShellCMD = insertionShellCMD + "oneinsert,chr1," + this_sv_begin + ",1:"
        elif diti_random <= diti_arr[2]:
            sv_type = "tandemRepeat"
            copy_num = random.randrange(2, 11)
            tandemRepeatShellCMD = tandemRepeatShellCMD + "onetandem_repeat,chr1," + this_sv_begin + "," + this_sv_end + "," + str(copy_num) + ",1:"
        else:
            sv_type = "inversion"
            inversionShellCMD = inversionShellCMD + "oneinversion,chr1," + this_sv_begin + "," + this_sv_end + ",1:"

    deletionShellCMD = deletionShellCMD[:-1]
    insertionShellCMD = insertionShellCMD[:-1]
    tandemRepeatShellCMD = tandemRepeatShellCMD[:-1]
    inversionShellCMD = inversionShellCMD[:-1]
    if deletionShellCMD == "":
        deletionShellCMD = "onedelete,chr1,100000000,100001000,1,0"
    if insertionShellCMD == "":
        insertionShellCMD = "oneinsert,chr1,100000000,1"
    if tandemRepeatShellCMD == "":
        tandemRepeatShellCMD = "onetandem_repeat,chr1,100000000,100001000,3,1"
    if inversionShellCMD == "":
        inversionShellCMD = "oneinversion,chr1,100000000,100001000,1"
    insert_len = random.randrange(5, 50)
    # snp= snv()
    command = 'python2 -B Wangshj.py -out /home/cloudam/ICGCsimulat_2/my_folder_'+ str(iturn) \
            + ' -readlen ' + str(this_sample_readLength) \
            + ' -depth ' + str(this_sample_depth) \
            + ' -delete ' + deletionShellCMD \
            + ' -tandem_repeat ' + tandemRepeatShellCMD \
            + ' -inversion ' + inversionShellCMD \
            # + ' -insert ' + insertionShellCMD \
            # + ' -fib ' + str(insert_len) \
            
    # if snp:  # 判断snp是否不为空
    #     command = command + ' -snv ' + str(snp)
    shFile.write(command.encode())
    #模拟结构变异结束


    shFile.close()

want_to_sim_number =100
variant_amount = str(50 + want_to_sim_number + want_to_sim_SV_number)
tsml_arr_convertINT = [int(10 * tsml_arr[0]), int(10 * (tsml_arr[1] - tsml_arr[0])), int(10 * (tsml_arr[2] - tsml_arr[1])), int(10 * (tsml_arr[3] - tsml_arr[2]))]
tsml_arr_convertINT_str = str(tsml_arr_convertINT[0]) + str(tsml_arr_convertINT[1]) + str(tsml_arr_convertINT[2]) + str(tsml_arr_convertINT[3])
diti_arr_convertINT = [int(10 * diti_arr[0]), int(10 * (diti_arr[1] - diti_arr[0])), int(10 * (diti_arr[2] - diti_arr[1])), int(10 * (diti_arr[3] - diti_arr[2]))]
diti_arr_convertINT_str = str(diti_arr_convertINT[0]) + str(diti_arr_convertINT[1]) + str(diti_arr_convertINT[2]) + str(diti_arr_convertINT[3])
# os.system("mv A_stableCallerPaperSimFlowShell.sh " + sampleID + "_" + this_sample_depth + "_" + this_sample_readLength + "_" + variant_amount + "_" + tsml_arr_convertINT_str + "_" + diti_arr_convertINT_str + "_" + str(SV_block) + ".sh")
# os.system("sh " + sampleID + "_" + this_sample_depth + "_" + this_sample_readLength + "_" + variant_amount + "_" + tsml_arr_convertINT_str + "_" + diti_arr_convertINT_str + "_" + str(SV_block) + ".sh")
# os.system("mv " + sampleID + "_" + this_sample_depth + "_" + this_sample_readLength + "_" + variant_amount + "_" + tsml_arr_convertINT_str + "_" + diti_arr_convertINT_str + "_" + str(SV_block) + ".sh A_batchSimuShells")
# os.system("mv " + sampleID + "_" + this_sample_depth + "_" + this_sample_readLength + " A_batchSimuSamples")
