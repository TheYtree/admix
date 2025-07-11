#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import numpy as np
from collections import defaultdict

def create_position_based_model(model_name):
    """创建基于位置的模型文件"""
    print(f"创建基于位置的模型: {model_name}")
    
    # 读取原始模型文件
    alleles_file = f"data/{model_name}.alleles"
    frequency_file = f"data/{model_name}.{admix_models.n_populations(model_name)}.F"
    
    if not os.path.exists(alleles_file) or not os.path.exists(frequency_file):
        print(f"错误: 模型文件不存在")
        return None
    
    # 创建位置到SNP的映射
    position_to_snp = {}
    
    # 读取alleles文件
    with open(alleles_file, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            if len(row) >= 3:
                rsid, minor_allele, major_allele = row[0], row[1], row[2]
                # 这里我们需要从外部数据库获取位置信息
                # 暂时使用rsID作为键
                position_to_snp[rsid] = {
                    'rsid': rsid,
                    'minor_allele': minor_allele,
                    'major_allele': major_allele
                }
    
    return position_to_snp

def create_position_based_reader():
    """创建基于位置的数据读取器"""
    
    def read_raw_data_by_position(data_file_name):
        """基于位置读取原始数据"""
        check_file(data_file_name)
        processed_data = {}
        position_data = {}
        
        with open(data_file_name, 'r') as data:
            data_reader = csv.reader(data, delimiter='\t')
            for row in data:
                if len(row) == 4:
                    rsid, chrom, pos, genotype = row
                    
                    # 检查基因型有效性
                    if genotype not in ['--', 'DD', 'II', 'DI'] and len(genotype) == 2:
                        if genotype[0] in ['A', 'T', 'G', 'C'] and genotype[1] in ['A', 'T', 'G', 'C']:
                            # 使用位置作为键
                            position_key = f"{chrom}:{pos}"
                            processed_data[position_key] = genotype
                            position_data[position_key] = {
                                'rsid': rsid,
                                'chrom': chrom,
                                'pos': pos,
                                'genotype': genotype
                            }
        
        return processed_data, position_data
    
    return read_raw_data_by_position

def create_position_based_model_reader():
    """创建基于位置的模型读取器"""
    
    def read_model_by_position(model):
        """基于位置读取模型数据"""
        # 这里我们需要创建一个包含位置信息的模型文件
        # 暂时返回空数据
        return [], [], [], []
    
    return read_model_by_position

# 修改现有的数据处理函数
def modify_raw_data_processing():
    """修改raw_data_processing.py以支持基于位置的匹配"""
    
    # 创建新的数据处理函数
    def twenty_three_and_me_by_position(data_file_name):
        check_file(data_file_name)
        processed_data = {}
        position_data = {}
        
        with open(data_file_name, 'r') as data:
            data_reader = csv.reader(data, delimiter='\t')
            for row in data:
                if len(row) == 4:
                    rsid, chrom, pos, genotype = row
                    
                    # 检查基因型有效性
                    if genotype not in ['--', 'DD', 'II', 'DI'] and len(genotype) == 2:
                        if genotype[0] in ['A', 'T', 'G', 'C'] and genotype[1] in ['A', 'T', 'G', 'C']:
                            # 使用位置作为键
                            position_key = f"{chrom}:{pos}"
                            processed_data[position_key] = genotype
                            position_data[position_key] = {
                                'rsid': rsid,
                                'chrom': chrom,
                                'pos': pos,
                                'genotype': genotype
                            }
        
        return processed_data
    
    return twenty_three_and_me_by_position

def create_position_based_admix_fraction():
    """创建基于位置的祖先成分计算函数"""
    
    def admix_fraction_by_position(model, raw_data_format, raw_data_file=None, tolerance=1e-3):
        """基于位置计算祖先成分"""
        # 读取基于位置的数据
        genome_data, position_data = read_raw_data_by_position(raw_data_format, raw_data_file)
        
        # 这里需要实现基于位置的SNP匹配和频率计算
        # 暂时返回均匀分布
        n_populations = admix_models.n_populations(model)
        return np.ones(n_populations) / n_populations
    
    return admix_fraction_by_position

if __name__ == "__main__":
    # 测试基于位置的读取
    test_file = "111-1975-6559.txt"
    
    # 创建基于位置的读取器
    read_raw_data_by_position = create_position_based_reader()
    
    # 测试读取
    try:
        processed_data, position_data = read_raw_data_by_position(test_file)
        print(f"成功读取基于位置的数据: {len(processed_data)} 个SNP")
        
        # 显示前几个位置
        count = 0
        for pos_key, data in position_data.items():
            if count < 5:
                print(f"{pos_key}: {data['genotype']}")
                count += 1
        
    except Exception as e:
        print(f"错误: {e}") 