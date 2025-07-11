#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv
import admix_models

def load_rsid_map(rsid_file):
    """加载rsid到(chromosome, position)的映射"""
    rsid_map = {}
    with open(rsid_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3 or row[0].startswith('#'):
                continue
            rsid, chrom, pos = row[0], row[1], row[2]
            rsid_map[rsid] = (chrom, pos)
    return rsid_map

def process_model(model, rsid_map, out_dir='data/position_models'):
    """为单个模型生成带位置信息的新模型文件，保持与原始文件行数一致"""
    alleles_file = f'data/{model}.alleles'
    out_file = os.path.join(out_dir, f'{model}.position.alleles')
    os.makedirs(out_dir, exist_ok=True)
    
    with open(alleles_file, 'r') as fin, open(out_file, 'w') as fout:
        reader = csv.reader(fin, delimiter=' ')
        writer = csv.writer(fout, delimiter='\t')
        
        for row in reader:
            if len(row) < 3:
                continue
            rsid, minor, major = row[0], row[1], row[2]
            
            if rsid in rsid_map:
                chrom, pos = rsid_map[rsid]
                writer.writerow([rsid, chrom, pos, minor, major])
            else:
                # 没有找到位置信息的rsid也写入，但chrom/pos留空
                writer.writerow([rsid, '', '', minor, major])
    
    print(f"模型 {model} 已生成: {out_file}")

def main():
    rsid_file = 'rsid.txt'
    rsid_map = load_rsid_map(rsid_file)
    models = admix_models.models()
    for model in models:
        process_model(model, rsid_map)

if __name__ == '__main__':
    main() 