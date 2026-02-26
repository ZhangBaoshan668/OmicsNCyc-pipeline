# OmicsNCyc：An integrated pipeline for fast analysis of nitrogen cycle functional genes with a curated database
Version：v1.0.0

Update：2026/02/26

## Background
To accurately analyze microbially driven nitrogen cycling, current databases and bioinformatics tools face critical limitations, as existing nitrogen-cycle databases (e.g., NCycDB) suffer from a lack of standardized analytical thresholds, leading to high false-positive rates, while also containing ambiguous taxonomic classifications and applying inappropriate uniform clustering cutoffs to functionally diverse genes such as amoA. Furthermore, although metagenomics overcomes the primer bias inherent in amplicon sequencing, available analysis pipelines remain fragmented, with current tools either designed for specific tasks or, like web-based platforms, lacking the flexibility required for large-scale or customized analyses. Therefore, a comprehensive solution is urgently needed, necessitating the development of a rigorously curated database with gene-specific thresholds for similarity and taxonomy, coupled with an integrated, user-friendly cloud platform and open-source pipeline to enable precise, high-throughput profiling of nitrogen-cycling communities.Therefore, we developed the OmicsNCyc pipeline.

## Pipeline manual and file description

Files description:

- README.md     # Introduction and install
- blast         # Directory for `BLAST` functional gene index files
- database      # Directory for `DIAMOND` functional gene index and reference database files
- software      # Directory for analysis software used in the pipeline
- sub           # Directory for all Python and R scripts used in the pipeline
- nitrogen_rpkm.py                   # Main script file

For Illumina or BGI-seq next-generation sequencing amplicon
- pipeline.sh # Command-line analysis for Windows and Linux
- pipeline_mac.sh # Command-line analysis for MacOS
- result/ # Example result data
- seq # short-read sequencing amplicon
- result/Diversity.Rmd # Interactive diversity analysis in R and output reproducible report in HTML format
- qiime2 # Using QIIME 2 analysis amplicon data 
- Accu16S_ITS # Absolute quantify amplicon analysis script
- advanced #

For PacBio or Nanopore third-generation long-read sequencing amplicon
- PacBio # Pipeline for PacBio long-read amplicon sequencing analysis
- Nanopore # Pipeline for Nanopore long-read amplicon sequencing analysis
- Mock # Synthetic community sequencing by Illumina and Pacbio for compare short and long amplicon
- snakemake # Pipeline for long-read sequencing amplicon





