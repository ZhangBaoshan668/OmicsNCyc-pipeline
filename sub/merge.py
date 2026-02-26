#!/usr/bin/python3
import os
import sys
import glob
import csv
import pandas as pd

def merge_csv_files(output_file, file_pattern, indir):
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        seen = set()
        csv_files = glob.glob(os.path.join(indir, '**', file_pattern), recursive=True)
        header_written = False  
        first_file = True  

        for csv_file in csv_files:
            with open(csv_file, 'r', newline='') as infile:
                reader = csv.reader(infile)
                header = next(reader, None)  
                if first_file:  
                    if header:  
                        writer.writerow(header)  
                        header_written = True
                for row in reader:
                    row_tuple = tuple(row)
                    if row_tuple not in seen:  
                        writer.writerow(row)  
                        seen.add(row_tuple)  
                first_file = False  

def convert_long_to_Count(input_file, output_file):
    # 读取CSV文件
    df = pd.read_csv(input_file)
    # 将长表格转换为宽表格
    wide_df = df.pivot(index='Taxon', columns='Sample ID', values='Count')
    wide_df.fillna(0, inplace=True)
    # 保存为新的CSV文件
    wide_df.to_csv(output_file)

def convert_long_to_wide(input_file, output_file):
    df = pd.read_csv(input_file)
    # 将长表格转换为宽表格
    wide_df = df.pivot(index='Taxon', columns='Sample ID', values='rpkm')
    wide_df.fillna(0, inplace=True)
    # 保存为新的CSV文件
    wide_df.to_csv(output_file)
# 获取命令行参数
indir = sys.argv[1]
gene  = sys.argv[2]
# 定义合并文件的输出目录
merge_dir = os.path.join(indir, 'merge')
if not os.path.exists(merge_dir):
    os.makedirs(merge_dir)


# 为特定基因创建子文件夹
gene_dir = os.path.join(merge_dir, gene)
if not os.path.exists(gene_dir):
    os.makedirs(gene_dir)

# 定义输出文件路径，确保它们位于merge目录中
result_output = os.path.join(gene_dir, f'{gene}_merged_results.txt')
target_output = os.path.join(gene_dir, f'{gene}_merged_target.fa')
rpkm_output = os.path.join(gene_dir, f'{gene}_merged_rpkm.txt')
tax_output = os.path.join(gene_dir, f'{gene}_merged_tax.txt')
total_reads_output = os.path.join(merge_dir, f'sample_merged_total_reads.txt')
kingdom_counts = os.path.join(gene_dir, f'{gene}_kingdom_merged_count.csv')
phylum_counts = os.path.join(gene_dir, f'{gene}_phylum_merged_count.csv')
class_counts = os.path.join(gene_dir, f'{gene}_class_merged_count.csv')
order_counts = os.path.join(gene_dir, f'{gene}_order_merged_count.csv')
family_counts = os.path.join(gene_dir, f'{gene}_family_merged_count.csv')
genus_counts = os.path.join(gene_dir, f'{gene}_genus_merged_count.csv')
species_counts = os.path.join(gene_dir, f'{gene}_species_merged_count.csv')


kingdom_cpm_output = os.path.join(gene_dir, f'{gene}_merged_kingdom_rpkm.csv')
phylum_cpm_output = os.path.join(gene_dir, f'{gene}_merged_phylum_rpkm.csv')
class_cpm_output = os.path.join(gene_dir, f'{gene}_merged_class_rpkm.csv')
order_cpm_output = os.path.join(gene_dir, f'{gene}_merged_order_rpkm.csv')
family_cpm_output = os.path.join(gene_dir, f'{gene}_merged_family_rpkm.csv')
genus_cpm_output = os.path.join(gene_dir, f'{gene}_merged_genus_rpkm.csv')
species_cpm_output = os.path.join(gene_dir, f'{gene}_merged_species_rpkm.csv')


# 合并 .result.txt 文件
with open(result_output, 'w') as outfile:
    txt_files = glob.glob(os.path.join(indir, '**', f'{gene}_*.result.txt'), recursive=True)
    for txt_file in txt_files:
        with open(txt_file, 'r') as infile:
            outfile.write(infile.read())

# 合并 .target.fa 文件
with open(target_output, 'w') as outfile:
    target_files = glob.glob(os.path.join(indir, '**', f'{gene}_*.target.fa'), recursive=True)
    for target_file in target_files:
        with open(target_file, 'r') as infile:
            outfile.write(infile.read())

# 合并 .rpkm.txt 文件
with open(rpkm_output, 'w') as outfile:
    rpkm_files = glob.glob(os.path.join(indir, '**', f'{gene}_*.rpkm.txt'), recursive=True)
    header_written = False
    # 遍历所有rpkm文件
    for rpkm_file in rpkm_files:
        with open(rpkm_file, 'r') as infile:
            # 逐行读取文件内容
            for line in infile:
                # 如果是第一个文件，则写入标题行
                if not header_written:
                    outfile.write(line)
                    header_written = True  # 标记标题行已写入
                else:
                    # 如果不是标题行，则写入数据行
                    if not line.startswith("SampleID"):  # 假设标题行以"Part"开头
                        outfile.write(line)
    
    
    
    # for rpkm_file in rpkm_files:
    #     with open(rpkm_file, 'r') as infile:
    #         outfile.write(infile.read())

# 合并 .tax.txt 文件
with open(tax_output, 'w') as outfile:
    tax_files = glob.glob(os.path.join(indir, '**', f'{gene}_*.tax_result.txt'), recursive=True)
    for tax_file in tax_files:
        with open(tax_file, 'r') as infile:
            outfile.write(infile.read())

# 合并 .total_reads.txt 文件
with open(total_reads_output, 'w') as outfile:
    total_files = glob.glob(os.path.join(indir, '**', '*_total_reads.txt'), recursive=True)
    for total_file in total_files:
        with open(total_file, 'r') as infile:
            outfile.write(infile.read())

# 合并分类层级的 .csv 文件
merge_csv_files(kingdom_counts, f'{gene}_*_kingdom_counts.csv', indir)
merge_csv_files(phylum_counts, f'{gene}_*_phylum_counts.csv', indir)
merge_csv_files(class_counts, f'{gene}_*_class_counts.csv', indir)
merge_csv_files(order_counts, f'{gene}_*_order_counts.csv', indir)
merge_csv_files(family_counts, f'{gene}_*_family_counts.csv', indir)
merge_csv_files(genus_counts, f'{gene}_*_genus_counts.csv', indir)
merge_csv_files(species_counts, f'{gene}_*_species_counts.csv', indir)

merge_csv_files(kingdom_cpm_output, f'{gene}_*_kingdom_rpkm.csv', indir)
merge_csv_files(phylum_cpm_output, f'{gene}_*_phylum_rpkm.csv', indir)
merge_csv_files(class_cpm_output, f'{gene}_*_class_rpkm.csv', indir)
merge_csv_files(order_cpm_output, f'{gene}_*_order_rpkm.csv', indir)
merge_csv_files(family_cpm_output, f'{gene}_*_family_rpkm.csv', indir)
merge_csv_files(genus_cpm_output, f'{gene}_*_genus_rpkm.csv', indir)
merge_csv_files(species_cpm_output, f'{gene}_*_species_rpkm.csv', indir)

convert_long_to_Count(kingdom_counts,kingdom_counts)
convert_long_to_Count(phylum_counts,phylum_counts)
convert_long_to_Count(class_counts,class_counts)
convert_long_to_Count(order_counts,order_counts)
convert_long_to_Count(family_counts,family_counts)
convert_long_to_Count(genus_counts,genus_counts)
convert_long_to_Count(species_counts,species_counts)


convert_long_to_wide(kingdom_cpm_output,kingdom_cpm_output)
convert_long_to_wide(phylum_cpm_output,phylum_cpm_output)
convert_long_to_wide(class_cpm_output,class_cpm_output)
convert_long_to_wide(order_cpm_output,order_cpm_output)
convert_long_to_wide(family_cpm_output,family_cpm_output)
convert_long_to_wide(genus_cpm_output,genus_cpm_output)
convert_long_to_wide(species_cpm_output,species_cpm_output)
