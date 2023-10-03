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
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1677814041&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Nzc4MTQwNDF9fX1dfQ__&Signature=CeQnToHDIRIiiKInKBCYNLXM3TnZOI346o1XXSjTiPCaAO-B4r1kVheGJF3ZkWWZG1ea1DJN9P4kQ9BHzrP--PHPYhatI5gIB8pxD74WoNAxa4iZxiQAYUKRq7I4z58L2jVwgTbddeOWfSSi2atH2zUVVOOAepnmGkH554a-zdBw2wT4CX1SABsFJ9uODCKkYn5wjg~KxOAK2ULT6FAWcf6IJWLM4OKm9Lz~ill576WFYIfl3OMng~tp-MPC-i641I9mb3qB4O1rsYSvI-IrvULbJ~aCUqNqW9pkdBZxfx5RHzv-SNpSau7WYLN613UMyKshMm07W-GFtwZ~XYYU5w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

```
To download dependencies that were developed internally by BWH Bioinformatics and Genomics Hub 
Note: Currently it is much easier to have conda handle the packages.

```
remotes::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")
remotes::install_github("bwh-bioinformatics-hub/qcreporter")
remotes::install_github('satijalab/seurat-wrappers')
git clone https://github.com/PAIN-initiative/rna_seq_pipeline_bwh.git
```
Github pages of dependencies developed internally: \
scRNA_seq Pipeline: https://github.com/PAIN-initiative/rna_seq_pipeline_bwh \
H5MANIPULATOR: https://github.com/PAIN-initiative/H5MANIPULATOR \
qcreporter: https://github.com/PAIN-initiative/qcreporter 

# 3. Nextflow installation 
To install nextflow:
```
git clone https://github.com/PAIN-initiative/nextflow-scRNAseq.git
cd nextflow-scRNAseq
cd nextflow
tar -xvzf nextflow-22.10.6.tar.gz
```

# 4. Conda 
To create conda environment with dependencies installed:
```
cd nextflow-scRNAseq
cd env/ 
conda env create -f environment.yml
conda activate scrna_nextflow_pipeline
Rscript install_R_packages.R
nextflow run scrna_seq.nf
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
