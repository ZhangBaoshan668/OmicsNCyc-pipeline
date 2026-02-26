# #!/usr/bin/python3
# import sys
# import os
# import pandas as pd
# from collections import defaultdict

# # 从命令行参数中获取输入
# gene = sys.argv[1]
# sample_id = sys.argv[2]
# blast_results_file = sys.argv[3]
# sample_dir = os.path.dirname(blast_results_file)  # 获取样本文件夹路径

# # 初始化分类层级的计数器
# taxonomic_counts = {
#     'Kingdom': defaultdict(int),
#     'Phylum': defaultdict(int),  
#     'Class': defaultdict(int),
#     'Order': defaultdict(int),
#     'Family': defaultdict(int),
#     'Genus': defaultdict(int),
#     'Species': defaultdict(int)
# }

# # 打开并逐行处理结果文件
# with open(blast_results_file, "r") as file:
#     for line in file:
#         # 跳过空行
#         if line.strip() == "":
#             continue
#         # 分割每一行的数据
#         parts = line.strip().split("\t")

#         # 检查是否包含第六列
#         if len(parts) < 6:
#             continue

#         # 获取第六列的相似性结果，并尝试转换为浮点数
#         similarity_str = parts[5].strip()
#         try:
#             similarity = float(similarity_str)
#         except ValueError:
#             continue

#         # 获取分类信息并分割成层级
#         taxonomic_info = parts[2]
#         taxonomic_levels = taxonomic_info.split(";")

#         # 初始化一个字典来存储当前行的分类层级
#         current_taxonomic_info = {
#             'Kingdom': None,
#             'Phylum': None,
#             'Class': None,
#             'Order': None,
#             'Family': None,
#             'Genus': None,
#             'Species': None
#         }

#         # 提取当前行的分类层级信息
#         for level in taxonomic_levels:
#             if level.startswith("d_"):  # Kingdom
#                 current_taxonomic_info['Kingdom'] = level[2:]
#             elif level.startswith("p_"):  # Phylum
#                 current_taxonomic_info['Phylum'] = level[2:]
#             elif level.startswith("c_"):  # Class
#                 current_taxonomic_info['Class'] = level[2:]
#             elif level.startswith("o_"):  # Order
#                 current_taxonomic_info['Order'] = level[2:]
#             elif level.startswith("f_"):  # Family
#                 current_taxonomic_info['Family'] = level[2:]
#             elif level.startswith("g_"):  # Genus
#                 current_taxonomic_info['Genus'] = level[2:]
#             elif level.startswith("s_"):  # Species
#                 current_taxonomic_info['Species'] = level[2:]

#         # 根据相似性调整分类信息
#         if gene == "AOA_amoA" and similarity < 81.48:
#             current_taxonomic_info['Genus'] = 'unclassified'
#             current_taxonomic_info['Species'] = 'unclassified'
#         elif gene == "AOA_amoA" and similarity < 94.91:
#             current_taxonomic_info['Species'] = 'unclassified'
#         elif gene == "AOA_amoB" and similarity < 59.06:
#             current_taxonomic_info['Genus'] = 'unclassified'
#             current_taxonomic_info['Species'] = 'unclassified'
#         elif gene == "AOA_amoB" and similarity < 92.40:
#             current_taxonomic_info['Species'] = 'unclassified'

#         # 统计每个分类层级的计数
#         for tax_level, taxon in current_taxonomic_info.items():
#             if taxon:  # 如果有分类信息，则统计
#                 taxonomic_counts[tax_level][taxon] += 1

# # 遍历每个分类层级的计数并保存到CSV文件
# for tax_level, counts in taxonomic_counts.items():
#     output_results = [[sample_id, taxon, count] for taxon, count in counts.items()]
#     df = pd.DataFrame(output_results, columns=['Sample ID', 'Taxon', 'Count'])
#     df.to_csv(os.path.join(sample_dir, f"{gene}_{sample_id}_{tax_level.lower()}_counts.csv"), index=False)

#!/usr/bin/python3
import sys
import os
import pandas as pd
from collections import defaultdict

# 定义基因特异性阈值 (放在文件开头便于维护)
GENE_THRESHOLDS = {
    "anfG":{"genus":50.43,"species":72.05},
    "ansB":{"genus":56.1,"species":82.32},
    "comammox_amoA":{"genus":100,"species":58.02},
    "comammox_nxrA":{"genus":48.83,"species":86.55},
    "AOA_amoA":{"genus":81.48,"species":94.91},
    "AOA_amoB":{"genus":59.06,"species":92.4},
    "AOA_amoC":{"genus":82.89,"species":82.94},
    "AOB_amoA":{"genus":42.49,"species":78.82},
    "AOB_amoB":{"genus":42.14,"species":75.95},
    "AOB_amoC":{"genus":54.8,"species":70.18},
    "hao":{"genus":53.58,"species":77.3},
    "hzsA":{"genus":62.58,"species":87.67},
    "nxrA":{"genus":23.12,"species":83.8},
    "nxrB":{"genus":31.09,"species":95.8},
    "glsA":{"genus":31.88,"species":44.17},
    "hcp":{"genus":15.23,"species":23.7},
    "hdh":{"genus":72.84,"species":76.43},
    "hox":{"genus":53.13,"species":64.36},
    "hzo":{"genus":74.45,"species":78.24},
    "hzsB":{"genus":62.97,"species":78.48},
    "hzsC":{"genus":74.93,"species":88.78},
    "napA":{"genus":52.29,"species":72.68},
    "napB":{"genus":32.21,"species":42.59},
    "napC":{"genus":43.93,"species":54.17},
    "narB":{"genus":28.15,"species":68.53},
    "narC":{"genus":34.96,"species":70.4},
    "narG":{"genus":48.21,"species":65.78},
    "narH":{"genus":51.36,"species":69.41},
    "narI":{"genus":31,"species":56.45},
    "narJ":{"genus":28.21,"species":42.29},
    "narV":{"genus":31,"species":56.43},
    "narW":{"genus":28.57,"species":42.86},
    "narY":{"genus":51.36,"species":69.41},
    "narZ":{"genus":47.11,"species":65.3},
    "nasA":{"genus":23.5,"species":62.45},
    "nasB":{"genus":37.92,"species":65.05},
    "nasD":{"genus":26.87,"species":47.18},
    "nasE":{"genus":33.43,"species":42.66},
    "nifD":{"genus":26.8,"species":27.07},
    "nifH":{"genus":44.88,"species":46.75},
    "nifK":{"genus":25.64,"species":26.72},
    "nifW":{"genus":31.3,"species":64.04},
    "nirA":{"genus":27.75,"species":44.94},
    "nirB":{"genus":36.29,"species":63.32},
    "nirD":{"genus":27.56,"species":46.73},
    "nirK":{"genus":22.4,"species":40.58},
    "nirS":{"genus":41.21,"species":56.49},
    "NOB_nxrA":{"genus":100,"species":83.61},
    "norB":{"genus":23.02,"species":43.62},
    "norC":{"genus":19.65,"species":62.16},
    "norD":{"genus":27.48,"species":33.61},
    "norQ":{"genus":25.86,"species":56.43},
    "nosZI":{"genus":47.5,"species":81.45},
    "nosZII":{"genus":35.12,"species":80.08},
    "nrfA":{"genus":30.02,"species":42.03},
    "nrfB":{"genus":30,"species":42.04},
    "nrfC":{"genus":38.43,"species":43.93},
    "nrfD":{"genus":28.12,"species":39.82},
    "nrfH":{"genus":27.94,"species":43.62},
    "pmoA":{"genus":44.66,"species":47.54},
    "pmoB":{"genus":36.1,"species":42.48},
    "pmoC":{"genus":43.57,"species":45.1},
    "ureA":{"genus":50,"species":52},
    "ureB":{"genus":48.54,"species":67.8},
    "ureC":{"genus":57.27,"species":64.14},
    "vnfD":{"genus":67.92,"species":98.52},
    "vnfG":{"genus":17.59,"species":89.8},
    "vnfH":{"genus":60.82,"species":93.45},
    "vnfK":{"genus":62.42,"species":93.05},
    "asnB":{"genus":27.10,"species":29.99},
    "glnA":{"genus":23.71,"species":25.71}
}

def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <gene> <sample_id> <blast_results_file>")
        sys.exit(1)

    gene = sys.argv[1]
    sample_id = sys.argv[2]
    blast_results_file = sys.argv[3]
    sample_dir = os.path.dirname(blast_results_file)

    # 初始化分类计数器
    taxonomic_counts = {
        'Kingdom': defaultdict(int),
        'Phylum': defaultdict(int),  
        'Class': defaultdict(int),
        'Order': defaultdict(int),
        'Family': defaultdict(int),
        'Genus': defaultdict(int),
        'Species': defaultdict(int)
    }

    with open(blast_results_file, "r") as file:
        for line in file:
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            try:
                similarity = float(parts[5].strip())
            except ValueError:
                continue

            # 解析分类信息
            current_taxonomic_info = parse_taxonomic_info(parts[2])
            
            # 应用基因特异性阈值
            apply_taxonomic_thresholds(gene, similarity, current_taxonomic_info)
            
            # 统计计数
            count_taxonomic_levels(current_taxonomic_info, taxonomic_counts)
    
    # 检查计数并应用 unclassified 条件
    apply_unclassified_condition(taxonomic_counts, 50)
    # 输出结果
    save_results(gene, sample_id, sample_dir, taxonomic_counts)

def parse_taxonomic_info(taxonomic_str):
    """解析分类信息字符串"""
    levels = {
        'd_': 'Kingdom',
        'p_': 'Phylum',
        'c_': 'Class',
        'o_': 'Order',
        'f_': 'Family',
        'g_': 'Genus',
        's_': 'Species'
    }
    
    current_info = {v: None for v in levels.values()}
    
    for level in taxonomic_str.split(";"):
        for prefix, tax_level in levels.items():
            if level.startswith(prefix):
                current_info[tax_level] = level[len(prefix):]
                break
                
    return current_info

def apply_taxonomic_thresholds(gene, similarity, tax_info):
    """应用基因特异性分类阈值"""
    thresholds = GENE_THRESHOLDS.get(gene, {})
    
    if not thresholds:
        return
        
    if similarity < thresholds.get("genus", float('inf')):
        tax_info['Genus'] = 'unclassified'
        tax_info['Species'] = 'unclassified'
    elif similarity < thresholds.get("species", float('inf')):
        tax_info['Species'] = 'unclassified'

def count_taxonomic_levels(tax_info, counters):
    """统计分类计数"""
    for level, taxon in tax_info.items():
        if taxon:
            counters[level][taxon] += 1

def apply_unclassified_condition(counters, threshold):
    """如果计数小于阈值，则标记为 unclassified"""
    for level, counts in counters.items():
        for taxon, count in list(counts.items()):
            if count < threshold:
                del counts[taxon]
                counts['unclassified'] += count
def save_results(gene, sample_id, output_dir, counters):
    """保存结果到CSV文件"""
    for tax_level, counts in counters.items():
        df = pd.DataFrame(
            [[sample_id, taxon, count] for taxon, count in counts.items()],
            columns=['Sample ID', 'Taxon', 'Count']
        )
        output_file = os.path.join(
            output_dir, 
            f"{gene}_{sample_id}_{tax_level.lower()}_counts.csv"
        )
        df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()