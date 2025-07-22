#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import numpy as np
import scipy.optimize as optimize
import sys
sys.path.append(os.path.dirname(__file__))
import admix_models

def check_file(file_name):
    if not os.path.isfile(file_name):
        print('Cannot find the file \'' + file_name + '\'!\n')
        exit()

def read_raw_data_by_position(data_file_name, data_format=None):
    check_file(data_file_name)
    processed_data = {}
    position_data = {}
    
    with open(data_file_name, 'r') as data:
        reader = csv.reader(data, delimiter='\t')
        
        # 读取前几行来检测数据格式
        sample_rows = []
        for i, row in enumerate(reader):
            if i < 5:  # 读取前5行来检测格式
                sample_rows.append(row)
            else:
                break
        
        # 检测数据格式
        detected_format = data_format
        if not detected_format:
            if len(sample_rows) > 0:
                if len(sample_rows[0]) == 5 and sample_rows[0][-1] in ['A', 'T', 'G', 'C']:
                    detected_format = 'ancestry'
                elif len(sample_rows[0]) == 4:
                    detected_format = '23andme'
                else:
                    detected_format = 'unknown'
        
        # 重新打开文件进行实际处理
        data.seek(0)
        reader = csv.reader(data, delimiter='\t')
        
        for row in reader:
            if detected_format == 'ancestry':
                # ancestry格式: rsid, chrom, pos, allele1, allele2
                if len(row) == 5 and row[-1] in ['A', 'T', 'G', 'C']:
                    rsid, chrom, pos, allele1, allele2 = row[0], row[1], row[2], row[3], row[4]
                    # 组合等位基因
                    genotype = allele1 + allele2
                    if len(genotype) == 2 and genotype[0] in ['A', 'T', 'G', 'C'] and genotype[1] in ['A', 'T', 'G', 'C']:
                        position_key = f"{chrom}:{pos}"
                        processed_data[position_key] = genotype
                        position_data[position_key] = {
                            'rsid': rsid,
                            'chrom': chrom,
                            'pos': pos,
                            'genotype': genotype
                        }
            else:
                # 默认23andme格式: rsid, chrom, pos, genotype
                if len(row) >= 4:
                    rsid, chrom, pos, genotype = row[0], row[1], row[2], row[3]
                    if genotype not in ['--', 'DD', 'II', 'DI'] and len(genotype) == 2:
                        if genotype[0] in ['A', 'T', 'G', 'C'] and genotype[1] in ['A', 'T', 'G', 'C']:
                            position_key = f"{chrom}:{pos}"
                            processed_data[position_key] = genotype
                            position_data[position_key] = {
                                'rsid': rsid,
                                'chrom': chrom,
                                'pos': pos,
                                'genotype': genotype
                            }
    
    return processed_data, position_data

def read_model_by_position(model_name):
    # 尝试多个可能的路径
    possible_paths = [
        f'data/position_models/{model_name}.position.alleles',
        f'admix/data/position_models/{model_name}.position.alleles',
        os.path.join(os.path.dirname(__file__), f'data/position_models/{model_name}.position.alleles')
    ]
    
    model_file = None
    for path in possible_paths:
        if os.path.isfile(path):
            model_file = path
            break
    
    if model_file is None:
        print(f'Cannot find the position model file for {model_name}!')
        print(f'Tried paths: {possible_paths}')
        exit()
    
    check_file(model_file)
    snp_keys = []
    minor_alleles = []
    major_alleles = []
    rsids = []
    with open(model_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 5:
                rsid, chrom, pos, minor, major = row[0], row[1], row[2], row[3], row[4]
                if chrom and pos:  # 只处理有位置信息的SNP
                    key = f"{chrom}:{pos}"
                    snp_keys.append(key)
                    minor_alleles.append(minor)
                    major_alleles.append(major)
                    rsids.append(rsid)
                else:  # 没有位置信息的SNP，添加空键
                    snp_keys.append('')
                    minor_alleles.append(minor)
                    major_alleles.append(major)
                    rsids.append(rsid)
    return snp_keys, minor_alleles, major_alleles, rsids

def genotype_matches_by_position(genome_data, snp_keys, minor_alleles, major_alleles):
    g_major = []
    g_minor = []
    matched_count = 0
    total_snps = len(snp_keys)
    
    for i, key in enumerate(snp_keys):
        if key and key in genome_data:  # 有位置信息且在数据中找到
            genotype = genome_data[key]
            major = major_alleles[i]
            minor = minor_alleles[i]
            major_count = sum(1 for allele in genotype if allele == major)
            minor_count = sum(1 for allele in genotype if allele == minor)
            g_major.append(major_count)
            g_minor.append(minor_count)
            matched_count += 1
        else:  # 没有位置信息或没找到匹配
            g_major.append(0)
            g_minor.append(0)
    
    print(f"位置匹配统计: {matched_count}/{total_snps} SNPs 匹配")
    return np.array(g_major), np.array(g_minor), matched_count, total_snps

def likelihood(g_major, g_minor, frequency, admixture_fraction):
    l1 = np.dot(g_major, np.log(np.dot(frequency, admixture_fraction)))
    l2 = np.dot(g_minor, np.log(np.dot(1 - frequency, admixture_fraction)))
    l = l1 + l2
    return -l

def admix_fraction_by_position(model, raw_data_file=None, tolerance=1e-3, data_format=None):
    genome_data, _ = read_raw_data_by_position(raw_data_file, data_format)
    print(f"读取到 {len(genome_data)} 个有效SNP")
    snp_keys, minor_alleles, major_alleles, rsids = read_model_by_position(model)
    
    # 尝试多个可能的频率文件路径
    freq_file_name = f"{model}.{admix_models.n_populations(model)}.F"
    possible_freq_paths = [
        f"data/{freq_file_name}",
        f"admix/data/{freq_file_name}",
        os.path.join(os.path.dirname(__file__), f"data/{freq_file_name}")
    ]
    
    freq_file = None
    for path in possible_freq_paths:
        if os.path.isfile(path):
            freq_file = path
            break
    
    if freq_file is None:
        print(f'Cannot find the frequency file for {model}!')
        print(f'Tried paths: {possible_freq_paths}')
        exit()
    
    frequency = []
    with open(freq_file, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            frequency.append([float(f) for f in row])
    frequency = np.array(frequency)
    g_major, g_minor, matched_count, total_snps = genotype_matches_by_position(genome_data, snp_keys, minor_alleles, major_alleles)
    n_populations = admix_models.n_populations(model)
    initial_guess = np.ones(n_populations) / n_populations
    bounds = tuple((0, 1) for i in range(n_populations))
    constraints = ({'type': 'eq', 'fun': lambda af: np.sum(af) - 1})
    likelihood_func = lambda af: likelihood(g_major, g_minor, frequency, af)
    admix_frac = optimize.minimize(
        likelihood_func,
        initial_guess,
        bounds=bounds,
        constraints=constraints,
        tol=tolerance).x
    
    # 计算利用率
    utilization_rate = matched_count / total_snps if total_snps > 0 else 0
    
    return admix_frac, utilization_rate

def print_results(model, result, sort=False, zh=False, ignore_zeros=False, utilization_rate=None):
    """打印分析结果，支持排序和中文显示"""
    populations = admix_models.populations(model)
    if not populations:
        print("无法获取种群信息")
        return
    
    admix_frac = np.array(result)
    populations = np.array(populations)
    
    # 执行降序排序
    if sort:
        idx = np.argsort(admix_frac)[::-1]
        admix_frac = admix_frac[idx]
        populations = populations[idx]
    
    # 在标题中显示利用率
    if utilization_rate is not None:
        print(f"\n=== {model} ({utilization_rate*100:.0f}%) 分析结果 ===")
    else:
        print(f"\n=== {model} 分析结果 ===")
    
    for i, (pop_en, pop_zh) in enumerate(populations):
        frac = admix_frac[i]
        if ignore_zeros and frac < 1e-4:
            continue
        if frac < 0.001:
            continue
        
        if zh:  # 中文显示
            population = f"{pop_zh} {pop_en}"
        else:  # 英文显示
            population = pop_en
        
        print(f"{population}: {frac*100:.2f}%")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='基于位置的祖先成分分析')
    parser.add_argument('-f', '--file', required=True, help='输入文件')
    parser.add_argument('-m', '--model', default='K7b', help='分析模型')
    parser.add_argument('-v', '--vendor', default='23andme', help='数据格式 (23andme, ancestry, ftdna, etc.)')
    parser.add_argument('-t', '--tolerance', default='1e-3', help='优化容差')
    parser.add_argument('--sort', action='store_true', help='对结果进行降序排序')
    parser.add_argument('-z', '--zhongwen', action='store_true', help='显示中文种群名称')
    parser.add_argument('--ignore-zeros', action='store_true', help='只显示非零比例')
    
    args = parser.parse_args()
    result_tuple = admix_fraction_by_position(args.model, args.file, float(args.tolerance), args.vendor)
    if isinstance(result_tuple, tuple) and len(result_tuple) == 2:
        result, utilization_rate = result_tuple
        if result is not None:
            print_results(args.model, result, args.sort, args.zhongwen, args.ignore_zeros, utilization_rate)
        else:
            print("分析失败")
    else:
        # 兼容旧版本返回格式
        if result_tuple is not None:
            print_results(args.model, result_tuple, args.sort, args.zhongwen, args.ignore_zeros)
        else:
            print("分析失败")

if __name__ == "__main__":
    main() 