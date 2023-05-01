#!/bin/bash
# Building the container requires the manual download of cellranger in the build directory
# From: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
sudo docker build .
