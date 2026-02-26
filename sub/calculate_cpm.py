#!/usr/bin/python3
import sys
import pandas as pd
import os

def calculate_cpm_rpkm(input_dir, sample_id, total_reads, gene_prefix, gene_length):
    # 定义分类层级的文件后缀
    taxonomic_suffixes = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    # 遍历每个分类层级
    for suffix in taxonomic_suffixes:
        input_file = os.path.join(input_dir, f"{gene_prefix}_{sample_id}_{suffix}_counts.csv")
        output_file = os.path.join(input_dir, f"{gene_prefix}_{sample_id}_{suffix}_rpkm.csv")
        
        # 检查输入文件是否存在
        if not os.path.exists(input_file):
            print(f"File not found: {input_file}")
            continue
        
        # 读取counts文件
        try:
            df = pd.read_csv(input_file)
        except pd.errors.EmptyDataError:
            print(f"No data in file: {input_file}")
            continue
        
        # 计算RPKM值并创建新的DataFrame
        rpkm_df = pd.DataFrame({
            'Sample ID': df['Sample ID'],  # 确保列名与CSV文件中的列名匹配
            'Taxon': df['Taxon'],  # 确保列名与CSV文件中的列名匹配
            'rpkm': ((df['Count'] / (total_reads*gene_length)) * 1e9).round(8)  # 除以基因长度
        })
        
        # 保存RPKM值到新的CSV文件
        rpkm_df.to_csv(output_file, index=False)
        # print(f"RPKM values for {suffix} have been calculated and saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 calculate_cpm.py <input_dir> <sample_id> <total_reads> <gene_prefix> <gene_length>")
        sys.exit(1)

    input_dir = sys.argv[1]
    sample_id = sys.argv[2]
    total_reads = int(sys.argv[3])
    gene_prefix = sys.argv[4]
    gene_length = float(sys.argv[5])  # 将基因长度转换为浮点数

    calculate_cpm_rpkm(input_dir, sample_id, total_reads, gene_prefix, gene_length)