cat gene.txt |cut -f1| xargs -I {} echo python3 /project/wujp/zbs/nitrogen_gene_pipline2/final_pipeline/sub/merge.py /project/wujp/zbs/Example_data {} >/project/wujp/zbs/Example_data/shell/mer_all.sh
sh /project/wujp/zbs/Example_data/shell/mer_all.sh
