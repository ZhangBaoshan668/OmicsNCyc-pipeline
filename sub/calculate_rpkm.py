# import sys
# import os
# def calculate_rpkm(total_reads_file, result_file, gene_file, output_dir, gene_name):
#     # 读取基因长度信息
#     gene_length = None
#     with open(gene_file, 'r') as f:
#         for line in f:
#             parts = line.strip().split()
#             if len(parts) == 2 and parts[0] == gene_name:
#                 gene_length = int(parts[1])
#                 break

#     if gene_length is None:
#         print(f"Length for gene {gene_name} not found in {gene_file}")
#         return

#     # 读取总读数
#     with open(total_reads_file, 'r') as total_f:
#         total_reads = int(total_f.readline().split()[1])
#         print(f"Total reads: {total_reads}")  # Debug print

#     # 提取结果文件中对应基因的读数并计算 RPKM 值
#     reads = 0
#     with open(result_file, 'r') as f:
#         for line in f:
#             reads = int(line.strip().split()[1])

#     if reads == 0:
#         print(f"Reads for gene {gene_name} not found in {result_file}")  # Debug print
#         return

#     rpkm = (reads / total_reads) * 1e9 / gene_length

#     # 构建输出文件路径，文件名为样本名.rpkm.txt
#     output_file_path = f"{output_dir}/{gene_name}.rpkm.txt"
#     # 将结果写入文件
#     with open(output_file_path, 'w') as out_f:   
#         out_f.write(f"{rpkm:.2f}\n")

# if __name__ == "__main__":
#     if len(sys.argv) != 6:
#         print("Usage: python calculate_rpkm.py <total_reads_file> <result_file> <gene_file> <output_dir> <gene_name>")
#         sys.exit(1)

#     total_reads_file = sys.argv[1]
#     result_file = sys.argv[2]
#     gene_file = sys.argv[3]
#     output_dir = sys.argv[4]
#     gene_name = sys.argv[5]  # 指定基因名称参数

#     # 确保输出目录存在
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     calculate_rpkm(total_reads_file, result_file, gene_file, output_dir, gene_name)
import sys
import os

def calculate_rpkm(total_reads_file, result_file, gene_file, output_dir, gene_name):
    # 读取基因长度信息
    gene_length = None
    with open(gene_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2 and parts[0] == gene_name:
                gene_length = int(parts[1])
                break

    if gene_length is None:
        print(f"Length for gene {gene_name} not found in {gene_file}")
        return

    # 读取总读数
    with open(total_reads_file, 'r') as total_f:
        total_reads = int(total_f.readline().split()[1])

    # 提取结果文件中对应基因的读数并计算 RPKM 值
    # reads = 0
    with open(result_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            reads = int(parts[1]) 

    rpkm = (reads / total_reads) * 1e9 / gene_length

    # 构建输出文件路径，文件名为样本名_基因名.rpkm.txt
    output_file_path = os.path.join(output_dir, f"{gene_name}_{parts[0]}.rpkm.txt")  # 使用样本名和基因名
    # 将结果写入文件
    with open(output_file_path, 'w') as out_f:
        out_f.write("SampleID\tRPKM\n")
        out_f.write(f"{parts[0]}\t{rpkm:.2f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python calculate_rpkm.py <total_reads_file> <result_file> <gene_file> <output_dir> <gene_name>")
        sys.exit(1)

    total_reads_file = sys.argv[1]
    result_file = sys.argv[2]
    gene_file = sys.argv[3]
    output_dir = sys.argv[4]
    gene_name = sys.argv[5]  # 指定基因名称参数

    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    calculate_rpkm(total_reads_file, result_file, gene_file, output_dir, gene_name)