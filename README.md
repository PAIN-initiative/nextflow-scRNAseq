# nextflow-scRNAseq
Nextflow for single-cell RNAseq

**************************
### Note: Currently please run using Conda Environment 
# 1. Introduction

### Overview:
Note: This pipeline is designed to be run post running CellRanger Count (cellrangers.nf). If the data meets basic QC from the CellRanger Websummary this pipeline will QC further and as well as run Seurat Preprocessing and Clustering.

This pipeline manages a scRNA-Seq workflow starting from raw fastq files and converting
them to standard file formats for use by downstream tools. The steps involved are:

* DoubletFinder: Single-Cell Remover of Doublets.
* Add meta, will add the metadata from CellRanger Count and user provided samplesheet to .h5 file.
* Create QC Report of all Samples provided.
* Cell Clustering and Cell Type Annotation.
* Perform Trajectory Analysis.
<a id="dependencies"></a>

# 2. Installation
## Dependencies    
This repository uses CellRanger Counts to generate the CellRangers outs directory that is used downstream: You can download it here:
CellRanger: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
```
You can Download CellRanger Software with this command: 
wget -O cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1697174417&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=U496cFRNxHuNHHpEKwFQgWNq8dUjTU~Vw3bLC35tlX9rLXeYNb1sCwMkP~GYXbsMFpt8MmbIYaFkpzGZoc1UuFNcH6-wL4wOpq-2EX6xkpo3wA3DSUoOf1vQzIvGXtRnLQCCKrQcei10~jO~WCS7xrmigrfP8sFTl0r7~AKzsNiMmPaiQflW~IyvocUm-r4lGeEim0twODUypxSneiT7ddhcDmrskb5v~NbNv9NTODMbfzgjAwVknPbVPgxPNhX9i3q-uf66IbNx7jQMVGRxkWjW2iG-6OglKO4M61RgNNL2UOEKdxUb18k5l7S0-J7gAOImUtp0xrZt1IiuX0HKTA__"


```
To download dependencies that were developed internally by BWH Bioinformatics and Genomics Hub 
Note: Currently it is much easier to have conda handle the packages.

``

# 3. Nextflow installation 
To install nextflow:
```
git clone https://github.com/PAIN-initiative/nextflow-scRNAseq.git
cd nextflow-scRNAseq
cd nextflow
tar -xvzf nextflow-22.10.6.tar.gz
cd nextflow-22.10.6/
sudo apt install openjdk-17-jre # if on linux
./nextflow

Nextflow will be successfully installed
```

# 4. Conda 
To create conda environment with dependencies installed:
```
cd nextflow-scRNAseq
cd env/ 
conda env create -f environment.yml
conda activate scrna_nextflow_pipeline
Rscript install_R_packages.R
pip install cellbender
cd ..
R # to activate R
install.packages('ggpubr')
# Leiden Algorithm Requirements
pip install leidenalg
pip install numpy
pip install pandas
nextflow/nextflow-22.10.6/nextflow run scrna_seq.nf
```

# 5. Containerization 

To deal with software dependencies and version controling a dockerfile has been created. \
To download docker image run 
```
docker pull acicalo4/snscrnaseq:latest
```

To mount data from local host to docker container run, example:
```
docker run -t -i -v path/to/data/you/want/mounted:/container/dir acicalo4/snscrnaseq /bin/bash
```

# 6. Example Usage: To Run After Pulling Docker Image
Setup:

Nextflow will parse a .csv file for the sample_ids and the path to the directory the fastq files are in for your project. Please provide at the minimum a sample_id column to the .csv file. \
If working with a .xls/.xlsx file please create a .csv file called samples.csv with a column labeled == 'sample_id' \
example: \

```
sample_id
P1708_SP093_105
```
```
nextflow run scrna_seq.nf -with-docker acicalo4/snscrnaseq
