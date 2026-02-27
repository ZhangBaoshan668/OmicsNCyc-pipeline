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
- `nitrogen_rpkm.py`                   # Main script file

## What can we do?

- Comprehensive analysis and visualization of N-cycling functional genes, including 71 genes such as *nifH*, *amoA* (bacterial and archaeal), *nirK*/*nirS*, *nosZ*, and *ureC*;
- From raw sequencing data to functional gene abundance tables and microbial community profiles;
- Accurate identification of N-cycling genes using gene-specific thresholds to minimize false positives;
- Analysis of a single N-cycling gene completes within five minutes.

![Figure 1](https://raw.githubusercontent.com/ABU789456/Figure/main/fig1.jpg)

**Figure 1. OmicsNCyc workflow. (A) Flowchart of major steps for NCycTaxDB construction. (B) Schematic of the data analysis workflow for nitrogen cycling functional genes. (C) Schematic overview of statistical and visual analyses via command-line or online web-platform modes of OmicsNCyc.**

## Main Features
+ Preprocessing and normalization of nitrogen cycling functional gene sequencing data
+ Functional gene abundance visualization
+ Taxonomic abundance visualization
+ Alpha diversity
+ Beta diversity
+ Differential abundance test
+ Network analysis

![Figure 2](https://raw.githubusercontent.com/ABU789456/Figure/main/fig2.jpg)

**Figure 2. Representative publication-ready visualizations.**

## Install

### Install Conda and R package

    # Create a new directory named “miniconda3” in your home directory.
    mkdir -p ~/miniconda3
    
    # Download the Linux Miniconda installation script for your chosen chip architecture and save the script as “miniconda.sh” in the miniconda3 directory.
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    
    # Run the “miniconda.sh” installation script in silent mode using bash.
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

    # Remove the “miniconda.sh” installation script file after installation is complete.
    rm ~/miniconda3/miniconda.sh

    # After installing, close and reopen your terminal application or refresh it by running the following command:
    source ~/miniconda3/bin/activate

    # Then, initialize conda on all available shells by running the following command:
    conda init --all

    # Add frequently used channels
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
    conda config --set show_channel_urls yes

    # Create Conda environment
    conda create -n OmicsNCyc
    conda activate OmicsNCyc

    # Conda install R
    conda install r-base=4.4.1

    # Install R package
    conda install r-BiocManager
    conda install r-ggplot2
    conda install r-devtools
    conda install r-vegan
    conda install r-igraph
    conda install r-dplyr
    conda install r-Hmisc
    conda install r-optparse
    conda install r-purrr
    conda install r-ggpubr
    conda install r-ggprism

    # Install amplicon package
    ## Enter R environment
    R
    library(devtools)
    install_github("microbiota/amplicon")
    ## Exit R environment
    q()

### Install OmicsNCyc

Install OmicsNCyc pipeline from github.

    # 

















