#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
sys.path.append(os.path.dirname(__file__))
from position_based_admix import admix_fraction_by_position, print_results
import admix_models

def batch_analysis(data_file, models, sort=False, zh=False, ignore_zeros=False, tolerance=1e-3):
    """批量分析多个模型"""
    print(f"开始批量分析，数据文件: {data_file}")
    print(f"分析模型: {', '.join(models)}")
    print("=" * 50)
    
    for model in models:
        try:
            result = admix_fraction_by_position(model, data_file, tolerance)
            if result is not None:
                print_results(model, result, sort, zh, ignore_zeros)
            else:
                print(f"\n=== {model} 分析失败 ===")
        except Exception as e:
            print(f"\n=== {model} 分析出错: {e} ===")
        
        print("-" * 30)

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='批量基于位置的祖先成分分析')
    parser.add_argument('-f', '--file', required=True, help='输入文件')
    parser.add_argument('-m', '--models', nargs='+', required=True, help='分析模型列表')
    parser.add_argument('-t', '--tolerance', default='1e-3', help='优化容差')
    parser.add_argument('--sort', action='store_true', help='对结果进行降序排序')
    parser.add_argument('-z', '--zhongwen', action='store_true', help='显示中文种群名称')
    parser.add_argument('--ignore-zeros', action='store_true', help='只显示非零比例')
    
    args = parser.parse_args()
    
    # 验证模型是否存在
    all_models = admix_models.models()
    for model in args.models:
        if model not in all_models:
            print(f"错误: 模型 {model} 不存在!")
            print(f"可用模型: {', '.join(all_models)}")
            return
    
    # 执行批量分析
    batch_analysis(
        args.file, 
        args.models, 
        args.sort, 
        args.zhongwen, 
        args.ignore_zeros, 
        float(args.tolerance)
    )

if __name__ == "__main__":
    main() 