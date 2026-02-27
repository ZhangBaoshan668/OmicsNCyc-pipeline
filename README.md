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
- gene length.txt    # Reference sequence length information of nitrogen cycle functional genes
- Example_data   # Example dataset

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

    # Create and activate OmicsNCyc environment
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

    # Create a new directory named “OmicsNCyc” in your home directory. 
    mkdir OmicsNCyc
    cd OmicsNCyc

    # Downdoald OmicsNCyc pipeline from github.
    git clone https://github.com/ABU789456/OmicsNCyc-pipeline.git

## Introduction to OmicsNCyc pipeline parameters

    # Parameter explanation
    cd OmicsNCyc-pipeline
    
    python3 nitrogen_rpkm.py -h
    usage: nitrogen_rpkm.py [-h] [-i CID] [-l LIS] [-g GENE] [-o OUTDIR] [-m GROUP] [-t THREADS] [-d ID] [-q QUERY] [-c CLUSTER] [-e EVALUE] [-v]
    
    options:
    -h, --help            show this help message and exit
    -i, --cid CID         input directory 
    -l, --lis LIS         sample list
    -g, --gene GENE       gene list
    -o, --outdir OUTDIR   output directory
    -m, --group GROUP     group list
    -t, --threads THREADS
                        number of threads to use for Parallel. Default is 10
    -d, --id ID           identity threshold for DIAMOND blastx. Default is 75
    -q, --query QUERY     query coverage threshold for DIAMOND blastx. Default is 75
    -c, --cluster CLUSTER
                        cluster threshold for USEARCH cluster_fast. Default is 0.97
    -e, --evalue EVALUE   E-value threshold for DIAMOND blastx. Default is 1e-5
    -v, --version         display version and author information and exit
**Note: The parameters -i, -l, -g, -m and -o are essential for the execution of the main script. The parameters -t, -d, -q, -c and -e are set as default values in the script, and their settings can also be adjusted according to specific needs.**

## Quick Start

    # Create work directory
    mkdir example
    cd example

    # Enter Conda environment
    conda activate OmicsNCyc

    # Create the "raw_data" folder in the working directory and the folder must be named "raw_data".
    # In the "raw_data" folder, only the single-end sequencing data of high-throughput sequencing should be stored, such as "_1.fastq.gz" or "_1.fq.gz".
    mkdir raw_data

    # Create "list.txt", "gene.txt", and "metadata.txt" file
    touch list.txt gene.txt metadata.txt
  
- **The `list.txt` file contains two columns of information. The first column represents the original sample number, and the second column represents the renamed sample number. The two columns are separated by a `tab`.**

![Figure 3](https://raw.githubusercontent.com/ABU789456/Figure/main/fig3.jpg)

- **The `metadata.txt` file must have a header row, with the columns labeled as `SampleID` and `Group`. The content of the header cannot be changed and it should contain two columns of information. The first column contains the renamed sample number, and the second column contains the sample grouping information. The two columns are separated by a `tab`.**

![Figure 4](https://raw.githubusercontent.com/ABU789456/Figure/main/fig4.jpg)

- **The `gene.txt` file contains two columns of information. The first column represents the gene name, and the second column represents the length of the gene reference sequence. The two columns are separated by a `tab`.**

![Figure 5](https://raw.githubusercontent.com/ABU789456/Figure/main/fig5.jpg)

**Note: The length information of the reference gene sequence is obtained from the `gene length.txt` file.**

**The `example` folder should contain the `raw_data` folder, `gene.txt`, `list.txt` and `metadata.txt`. Then, by running the main script file `nitrogen_rpkm.py` in the `example` folder, the analysis can be started.**

![Figure 6](https://raw.githubusercontent.com/ABU789456/Figure/main/fig6.jpg)

    # running
    
    python3 /your_path/nitrogen_rpkm.py -i Example_data -l list.txt -g gene.txt -m metadata.txt -o ./





















