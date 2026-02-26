#!/usr/bin/python3

import argparse
import os
import glob
import re
import threading

# Define version number and author information
VERSION = "1.0.0"
AUTHOR = "MEG"
LOGO = r'''
******************************************************************************                  
                             .%%%%%%%%%%%%%%%%%%%%%                             
                        %%%%% %%%%,.%% %%%  % %%%%% %%%%%                       
                     %%%%%%% % % %%%% %%%%%%%%%% % %%# /%%%%                    
                  %%% %% %.%%%%%%#             %%%%%#,% .%%%%%%                 
                %%%%  % %%%%                         %%%%,%%%%%%%               
              %%# #%%%%%%         (((((( ((((((         %%%% %%%%%%             
            %%%%#  %%%       #####(#(*     /((((((((       %%%%%  %%%           
           %%%%  %%%%     ######                 (((((      %%%  %%%%%          
          %%% %%%%%        (((((/                   ((##      %% % %%%%         
         #%%%  %%%     ((((((   ////(         ((((    ####     %% %%%%%         
         %%%%%%%%%    ((((       (((((((      ((((     ####    %%%  %%%%        
         %%%%%.%%     (((        ((((((((     ((#(     ####     %% %%%%%        
         %%%  %%%     #((        ((((  ((((   ####      ###     %%(  %%%        
         %%%%%#%%     ###        ((((   ((#(  ####     %##%     %%.% %%%        
         %%  %%%%%    ####       (##(    ,########     %%%     %%%% %%%%        
         ,%% %  %%     ####      ####      #######    %%%%     %%%%%%%%         
          %%%%%%%%%      ####    ####,       #####   %%       %%%%%%%%%         
           %%%%%%%%%%      #####                  #%%%      %%%%%%%%%%          
            %%%%%%%%%%*      #########/   #%%%%%%%%%      %%%%%%%%%%#           
              %%%%%%%%%%%         #%%%%%%%%%%%%         %%%%%%%%%%#             
                %%%%*%%%%%%%                         %%%%%%*%*%%%               
                  %%%%%%**%%%%%%%%             %%%%%%%(**%#%%%%                 
                     %%%%%***%*%%%%*%%%%%%%%#*%%%**#**%%%%%%                    
                        *%%%%%(%%**%*%**%*%*%%%**%%%%%%%                        
                              %%%%%%%%%%%%%%%%%%%%%         
Image design: Nancy Merino (2018);
ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/                                           
******************************************************************************
                              '''
parser = argparse.ArgumentParser(
    description=LOGO,
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('-i', '--cid', required=False, type=str, help='input directory')
parser.add_argument('-l', '--lis', required=False, type=str, help='sample list')
parser.add_argument('-g', '--gene', required=False, type=str, help='gene list')
parser.add_argument('-o', '--outdir', required=False, type=str, help='output directory')
parser.add_argument('-m', '--group', required=False, type=str, help='group list')
parser.add_argument('-t', '--threads', type=int, default=36, help='number of threads to use for Parallel. Default is 10')
parser.add_argument('-d', '--id', type=float, default=75, help='identity threshold for DIAMOND blastx. Default is 75')
parser.add_argument('-q', '--query', type=float, default=75, help='query coverage threshold for DIAMOND blastx. Default is 75')
parser.add_argument('-c', '--cluster', type=float, default=0.97, help='cluster threshold for USEARCH cluster_fast. Default is 0.97')
parser.add_argument('-e', '--evalue', default=1e-5, help='E-value threshold for DIAMOND blastx. Default is 1e-5')
parser.add_argument('-v', '--version', action='store_true', help='display version and author information and exit')

args = parser.parse_args()
if args.version:
    print(f"Version: {VERSION}")
    print(f"Author: {AUTHOR}")
    exit(0)

chip_id = args.cid
sample_list = args.lis
gene_list = args.gene
outdir = args.outdir
threads = args.threads
identity = args.id
query_cover = args.query
evalue = args.evalue
group_list = args.group
cluster = args.cluster


sample_list = os.path.abspath(sample_list)
outdir = os.path.abspath(outdir)
out_path = outdir
if not os.path.exists(out_path):
    os.mkdir(out_path)

script_dir = os.path.dirname(os.path.realpath(__file__))
sub_dir = os.path.join(script_dir, "sub")
db_dir = os.path.join(script_dir, "database")
blast_dir = os.path.join(script_dir, "blast")
software_dir = os.path.join(script_dir, "software")

# 读取文件重命名列表
rename_dict = {}
if os.path.exists(sample_list):
    with open(sample_list, 'r') as file:
        for line in file:
            original, new = line.strip().split()
            rename_dict[original] = new

fq_path = os.path.join(out_path, "raw_data")
if not os.path.exists(fq_path):
    os.mkdir(fq_path)

shell_path = os.path.join(out_path, "shell")
if not os.path.exists(shell_path):
    os.mkdir(shell_path)

all_samp = []
for index in rename_dict: 
    samp = rename_dict[index]
    all_samp.append(samp)
    samp_path = os.path.join(out_path, samp)
    if not os.path.exists(samp_path):
        os.mkdir(samp_path)
    possible_paths = glob.glob(os.path.join(fq_path, samp + '*.gz'))
    # 检查是否有匹配的文件
    if possible_paths:
        # 获取第一个匹配的文件路径
        full_samp_path = possible_paths[0]
        # 重命名文件
        new_full_samp_path = possible_paths[0].replace(samp, rename_dict[index])
        os.rename(full_samp_path, new_full_samp_path)
    else:
        print(f"No matching file found for {samp} in {fq_path}")
        continue
    with open(gene_list, 'r') as file:
        next(file)
        for line in file:
            gene, length = line.strip().split()
    with open(os.path.join(shell_path, samp+".sh"), 'w') as sh_file:
        sh_file.write(f'{software_dir}/fastp -i {new_full_samp_path} -o {samp_path}/{samp}.clean.fq.gz -w 16 -j {samp_path}/{samp}.json\n')
        sh_file.write(f'echo -n "{samp}" | paste - - <(jq \'.summary.after_filtering.total_reads\' {samp_path}/{samp}.json) | awk \'{{print $1 "\t" $2}}\' > {samp_path}/{samp}_total_reads.txt\n')
        sh_file.write(f'{software_dir}/seqtk seq -a  {samp_path}/{samp}.clean.fq.gz  >{samp_path}/{samp}.fa\n')
        sh_file.write(f'awk \'BEGIN{{FS=""; OFS=""; id=1}} /^>/ {{print ">{samp}_"id++; next}} {{print}}\' {samp_path}/{samp}.fa > {samp_path}/{samp}.rename.fa\n')
        sh_file.write(f'cat {gene_list} |cut -f1|xargs -I {{}} echo {software_dir}/diamond blastx -d {db_dir}/{{}} -q {samp_path}/{samp}.rename.fa -o {samp_path}/{{}}_{samp}.txt -e {evalue} --id {identity} --query-cover {query_cover} -f 6 -p 140 -k 1 > {samp_path}/diamond_blas.sh\n')    
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/diamond_blas.sh {threads}\n')
        sh_file.write(f"""cat {gene_list} |cut -f1|xargs -I {{}} echo 'paste -d '"'\\t'"' <(echo -n "{samp}") <(cat {samp_path}/{{}}_{samp}.txt | cut -f1 | sort | uniq | wc -l) > {samp_path}/{{}}_{samp}.result.txt' > {samp_path}/sort_uniq.sh\n""")
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/sort_uniq.sh {threads}\n')
        sh_file.write(f'cat {gene_list} |cut -f1|xargs -I {{}} echo "cat {samp_path}/{{}}_{samp}.txt | cut -f1 | sort | uniq > {samp_path}/{{}}_{samp}.id" > {samp_path}/id.sh\n')
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/id.sh {threads}\n')
        sh_file.write(f'cat {gene_list} |cut -f1|xargs -I {{}} echo "{software_dir}/seqtk subseq {samp_path}/{samp}.rename.fa {samp_path}/{{}}_{samp}.id > {samp_path}/{{}}_{samp}.target.fa" >{samp_path}/target.sh\n')
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/target.sh {threads}\n')
        sh_file.write(f'cat {gene_list} |cut -f1|xargs -I {{}} echo python3 {sub_dir}/calculate_rpkm.py {samp_path}/{samp}_total_reads.txt {samp_path}/{{}}_{samp}.result.txt {gene_list} {samp_path} {{}} >{samp_path}/rpkm.sh\n')
        sh_file.write(f'sh {samp_path}/rpkm.sh\n')
        sh_file.write(f'cat {gene_list} |cut -f1| xargs -I {{}} echo {software_dir}/blastx -db {blast_dir}/{{}} -query {samp_path}/{{}}_{samp}.target.fa -out {samp_path}/{{}}_{samp}.tax_result.txt -max_target_seqs 1 -outfmt \\\'6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname\\\' -num_threads 50 -evalue 1e-3 > {samp_path}/blastx.sh\n')
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/blastx.sh {threads}\n')
        sh_file.write(f'cat {gene_list} |cut -f1| xargs -I {{}} echo python3 {sub_dir}/tax_statistics.py {{}} {samp} {samp_path}/{{}}_{samp}.tax_result.txt > {samp_path}/tax_result.sh\n')
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/tax_result.sh {threads}\n')
        sh_file.write(f"""cat {gene_list} |cut -f1| xargs -I {{}} echo "cat {gene_list} |cut -f2|xargs -I {{[]}} sh -c \\"cat {samp_path}/{samp}_total_reads.txt|cut -f2|xargs -I {[]} python3 {sub_dir}/calculate_cpm.py {samp_path} {samp} {[]} {{}} {{[]}}\\"" > {samp_path}/cpm_commands.sh\n""")
        sh_file.write(f'python3 {sub_dir}/ParallelShellExecutor.py {samp_path}/cpm_commands.sh {threads}\n')
        sh_file.write(f'rm -rf {samp_path}/{samp}.clean.fq.gz\n')
        sh_file.write(f'rm -rf {samp_path}/{samp}.fa\n')
        sh_file.write(f'rm -rf {samp_path}/{samp}.rename.fa\n')
        sh_file.write(f'rm -rf {samp_path}/*plog\n')
        sh_file.write(f'rm -rf {samp_path}/*sh\n')

        
with open(os.path.join(shell_path, "all.sample.sh"), 'w') as all_sh:
    for samp in all_samp:
        all_sh.write(f'sh {shell_path}/{samp}.sh\n')

os.system(f'python3 {sub_dir}/ParallelShellExecutor.py {shell_path}/all.sample.sh {threads}') 

with open(os.path.join(shell_path, "merge.sh"), 'w') as sh_merge:
    sh_merge.write(f'cat {gene_list} |cut -f1| xargs -I {{}} echo python3 {sub_dir}/merge.py {out_path} {{}} >{shell_path}/mer_all.sh\n')
    sh_merge.write(f'sh {shell_path}/mer_all.sh\n')

os.system(f'sh {shell_path}/merge.sh')

with open(os.path.join(shell_path, "feature.sh"), 'w') as otutab:
    otutab.write(f"""cat {gene_list} | cut -f1 | xargs -I {{}} sh -c 'sed "s/_[0-9][0-9]*$//" {out_path}/merge/{{}}/{{}}_merged_target.fa > {out_path}/merge/{{}}/{{}}_merged_target_rename.fa'\n""")
    otutab.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -fastx_uniques {out_path}/merge/{{}}/{{}}_merged_target_rename.fa -fastaout {out_path}/merge/{{}}/{{}}_uniques.fa -sizeout -relabel OTU_  -minuniquesize 2\n')
    otutab.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -cluster_fast {out_path}/merge/{{}}/{{}}_uniques.fa -centroids {out_path}/merge/{{}}/{{}}_otus.fa -uc {out_path}/merge/{{}}/{{}}_clusters.uc -id {cluster} -minsize 2\n')
    otutab.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/vsearch --usearch_global {out_path}/merge/{{}}/{{}}_merged_target_rename.fa --db {out_path}/merge/{{}}/{{}}_otus.fa --id  {cluster} --threads 140  --otutabout {out_path}/merge/{{}}/{{}}_otutab.txt\n')
os.system(f'sh {shell_path}/feature.sh') 

with open(os.path.join(shell_path, "diversity.sh"), 'w') as diversity:
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} mkdir -p {out_path}/merge/{{}}/alpha {out_path}/merge/{{}}/beta {out_path}/merge/{{}}/network {out_path}/merge/{{}}/taxonomy\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/otutab_rare.R --input {out_path}/merge/{{}}/{{}}_otutab.txt --normalize {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.txt --output {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -otutab_stats {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.txt -output {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.stat\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -alpha_div_rare {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.txt -output {out_path}/merge/{{}}/alpha/{{}}_alpha_rare.txt -method without_replacement\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} sed -i "s/-/\t0.0/g" {out_path}/merge/{{}}/alpha/{{}}_alpha_rare.txt\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -beta_div {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.txt  -filename_prefix {out_path}/merge/{{}}/beta/\n')
    diversity.write(f'cat {gene_list} |cut -f1| xargs -I {{}} {software_dir}/usearch -otutab_counts2freqs {out_path}/merge/{{}}/alpha/{{}}_otutab_rare.txt -output {out_path}/merge/{{}}/{{}}_otutab_freqs.txt\n')
os.system(f'sh {shell_path}/diversity.sh')

with open(os.path.join(shell_path, "plot.sh"), 'w') as plot:
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_rare_curve.R --input {out_path}/merge/{{}}/alpha/{{}}_alpha_rare.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index richness --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index chao1 --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index shannon --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index simpson --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index invsimpson --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/alpha_boxplot.R --alpha_index ACE --input {out_path}/merge/{{}}/alpha/{{}}_alpha_div.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/alpha/ --width 89 --height 59\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/beta_pcoa.R --input {out_path}/merge/{{}}/beta/bray_curtis.txt --design {group_list} --group Group --output {out_path}/merge/{{}}/beta/bray_curtis.pcoa.pdf --width 89 --height 59 --label FALSE\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_phylum_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_phylum_rpkm --width 89 --height 95\n')    
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_class_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_class_rpkm --width 89 --height 95\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_order_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_order_rpkm --width 89 --height 95\n')    
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_family_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_family_rpkm --width 89 --height 95\n')    
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_genus_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_genus_rpkm --width 89 --height 95\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/tax_stackplot.R --input {out_path}/merge/{{}}/{{}}_merged_species_rpkm.csv --design {group_list} --group Group --color Paired --legend 9 --output {out_path}/merge/{{}}/taxonomy/{{}}_species_rpkm --width 89 --height 95\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/gene_rpkm_boxplot.R --input {out_path}/merge/{{}}/{{}}_merged_rpkm.txt --group {group_list}  --output {out_path}/merge/{{}}/{{}}_rpkm.pdf --width 160 --height 110\n')
    plot.write(f'cat {gene_list} |cut -f1| xargs -I {{}} Rscript {sub_dir}/network_plot.R --input {out_path}/merge/{{}}/{{}}_otutab.txt --group {group_list}  --output {out_path}/merge/{{}}/network --width 16 --height 11\n')

os.system(f'sh {shell_path}/plot.sh')

print ("Pipeline finished")



