# nextflow-scRNAseq
Nextflow for single-cell RNAseq

**************************

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
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1677814041&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Nzc4MTQwNDF9fX1dfQ__&Signature=CeQnToHDIRIiiKInKBCYNLXM3TnZOI346o1XXSjTiPCaAO-B4r1kVheGJF3ZkWWZG1ea1DJN9P4kQ9BHzrP--PHPYhatI5gIB8pxD74WoNAxa4iZxiQAYUKRq7I4z58L2jVwgTbddeOWfSSi2atH2zUVVOOAepnmGkH554a-zdBw2wT4CX1SABsFJ9uODCKkYn5wjg~KxOAK2ULT6FAWcf6IJWLM4OKm9Lz~ill576WFYIfl3OMng~tp-MPC-i641I9mb3qB4O1rsYSvI-IrvULbJ~aCUqNqW9pkdBZxfx5RHzv-SNpSau7WYLN613UMyKshMm07W-GFtwZ~XYYU5w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

This repository requires that `pandoc` and `libhdf5-devel` libraries are installed as dependencies of the `H5MANIPULATOR` functions:
```
sudo apt-get install pandoc libhdf5-dev
```

CRAN packages can be installed in R using:
```
install.packages("pacman",repo="https://cloud.r-project.org")
install.packages("remotes",repo="https://cloud.r-project.org")
install.packages('viridis',repo="https://cloud.r-project.org")
install.packages('cowplot',repo="https://cloud.r-project.org")
install.packages("scCustomize",repo="https://cloud.r-project.org")
install.packages("gt",repo="https://cloud.r-project.org")
install.packages("rio",repo="https://cloud.r-project.org")

# Dot plot is depedent on GitHub Report (https://github.com/Simon-Leonard/FlexDotPlot)
remotes::install_github("Simon-Leonard/FlexDotPlot")
```
Some Packages are Dependent on BiocManager
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scry")
BiocManager::install("ComplexHeatmap")
BiocManager::install("rhdf5")

```
To download dependencies that were developed internally by BWH Bioinformatics and Genomics Hub 
Note: It is much easier to pull the docker image.

```
remotes::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")
remotes::install_github("bwh-bioinformatics-hub/qcreporter")
remotes::install_github('satijalab/seurat-wrappers')
git clone https://github.com/bwh-bioinformatics-hub/rna_seq_pipeline_bwh.git
```
Github pages of dependencies developed internally: \
scRNA_seq Pipeline: https://github.com/PAIN-initiative/rna_seq_pipeline_bwh \
H5MANIPULATOR: https://github.com/PAIN-initiative/H5MANIPULATOR \
qcreporter: https://github.com/PAIN-initiative/qcreporter 

To create conda environment with dependencies install
```
conda env create -f environment.yml 
```

# 3. Containerization 

To deal with software dependencies and version controling a dockerfile has been created. \
To download docker image run 
```
docker pull acicalo4/snscrnaseq:latest
```

To mount data from local host to docker container run, example:
```
docker run -t -i -v path/to/data/you/want/mounted:/container/dir acicalo4/snscrnaseq /bin/bash
```

# 4. Example Usage: To Run After Pulling Docker Image
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
