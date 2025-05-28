
import argparse
import os

import re

def extract_variant_info(vcf_path):
    variant_data = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            # Skip metadata and header lines
            if line.startswith("#"):
                continue
            
            fields = line.strip().split('\t')
            
            # Ensure the line has enough fields
            if len(fields) < 8:
                continue
            
            # Chromosome and position
            chrom = fields[0]
            pos = fields[1]
            
            # Variant type (from ALT column)
            variant_type = fields[4]
            
            # Extract END from INFO field
            info_fields = fields[7].split(';')
            end_position = None
            for info in info_fields:
                if info.startswith("END="):
                    end_position = info.split('=')[1]
                    break
            
            # If END is not found in INFO, try to extract from ALT
            if end_position is None:
                alt = fields[4]
                match = re.search(r'\[(\d+):(\d+)\[|\](\d+):(\d+)\]', alt)
                if match:
                    end_position = match.group(2) or match.group(4)
            
            variant_data.append((chrom, pos, end_position, variant_type))
    
    return variant_data

def turth_info(vcf_path):
    variant_data = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            # Skip metadata and header lines
            if line.startswith("#"):
                continue
            
            fields = line.strip().split('\t')
            
            # Ensure the line has enough fields
            if len(fields) < 8:
                continue
            
            # Chromosome and position
            chrom = fields[0]
            pos = fields[1]
            
            # Variant type (from ALT column or INFO field)
            alt = fields[4]
            info_fields = fields[7].split(';')
            
            # Default variant type
            variant_type = "UNK"

            # 判断是否是结构变异
            if alt.startswith("<") and alt.endswith(">"):
                variant_type = alt  # Extract type from ALT field, e.g., <DEL>, <INV>
            else:
                for info in info_fields:
                    if info.startswith("SVTYPE="):
                        variant_type = info.split('=')[1]
                        break

            # Extract END from INFO field
            end_position = pos  # Default end position is the start position
            for info in info_fields:
                if info.startswith("END="):
                    end_position = info.split('=')[1]
                    break
            
            variant_data.append((chrom, pos, end_position, variant_type))
    
    return variant_data

def sort_vcf(input_path, output_path):
    # Extract header, metadata, and variants
    headers, variants = [], []
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line.strip())
            elif line.startswith("#"):
                column_headers = line.strip()
            else:
                variants.append(line.strip())

    # Sort variants by chromosome and then by position
    sorted_variants = sorted(variants, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))

    # Write to output
    with open(output_path, 'w') as f:
        for header in headers:
            f.write(header + "\n")
        f.write(column_headers + "\n")
        for variant in sorted_variants:
            f.write(variant + "\n")

def is_matching_variant(var1, var2):
    chrom1, start1, end1, type1 = var1
    chrom2, start2, end2, type2 = var2

    # Convert to integer if not None, otherwise keep as None
    start1 = int(start1) if start1 is not None else None
    end1 = int(end1) if end1 is not None else None
    start2 = int(start2) if start2 is not None else None
    end2 = int(end2) if end2 is not None else None

    if chrom1 == chrom2:
        start_condition = 0.25 * (1 - abs(start1 - start2) / max(start1, start2)) if start1 is not None and start2 is not None else 0
        end_condition = 0.25 * (1 - abs(end1 - end2) / max(end1, end2)) if end1 is not None and end2 is not None else 0
        score = 0.5 + start_condition + end_condition
    else:
        score = 0
    
    return score

def compute_f1(standard_vcf, called_vcf, threshold=0.5):
    TP = FN = FP = 0
    matched_indices = set()
    for std_idx, std_var in enumerate(standard_vcf):
        best_score = 0
        best_match_idx = None
        for call_idx, call_var in enumerate(called_vcf):
            score = is_matching_variant(std_var, call_var)
            if score > best_score:
                best_score = score
                best_match_idx = call_idx

        if best_score >= threshold:
            matched_indices.add(best_match_idx)
            TP += best_score
            FN += 1-best_score
        else:
            FN += 1

    FP = len(called_vcf) - len(matched_indices)

    precision = TP / (TP + FP) if (TP + FP) != 0 else 0
    recall = TP / (TP + FN) if (TP + FN) != 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) != 0 else 0

    return TP, FP, FN, f1, precision, recall

def evamain(input_vcf_1, input_vcf_2):
    # 排序输入VCF文件
    # sorted_vcf_1 = input_vcf_1 + ".sorted"
    # sorted_vcf_2 = input_vcf_2 + ".sorted"
    # sort_vcf(input_vcf_1, sorted_vcf_1)
    # sort_vcf(input_vcf_2, sorted_vcf_2)
    
    # 提取变异信息
    variant_info_1 = extract_variant_info(input_vcf_1)
    variant_info_2 = extract_variant_info(input_vcf_2)
    
    # 计算F1分数
    TP, FP, FN, f1_score, precision, recall = compute_f1(variant_info_1, variant_info_2)

    return TP, FP, FN, f1_score, precision, recall

