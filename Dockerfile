FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ARG R_VERSION=4.2.3
ARG OS_IDENTIFIER=ubuntu-2204

## Resolving R and lib dependencies
RUN apt-get update \
    && apt-get install -y \
    vim \
    wget \
    default-jdk \
    git \ 
    build-essential \
    python3 \
    python3-pip \
    libbz2-dev \
    libc6-dev \
    libgcc-9-dev \
    gcc-9-base \
    liblzma-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    pandoc \
    libpng-dev \
    pkg-config 

# Python Dependencies 
RUN pip install --no-cache-dir --upgrade pip && \
    pip install scrublet
# Install R
RUN wget https://cdn.posit.co/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*

## Install packages from CRAN & BioConductor
RUN Rscript -e "install.packages(c('BiocManager', 'devtools','FactoMineR', 'ggforce','systemfonts', 'assertthat', 'cowplot', 'data.table', 'dplyr', 'ids', 'ggplot2', 'jsonlite', 'Matrix', 'optparse', 'GetoptLong', 'purrr','R.utils','rmarkdown','stringr','gt','plotly','tidyr','Seurat','future','future.apply','rio','purrr','scCustomize','egg','DT','SoupX','reticulate','glmpca','remotes','viridis','qs','gridExtra','plyr','circlize','naniar','ggpubr','XML','RCurl'),repos='http://cran.us.r-project.org'); BiocManager::install(c('rhdf5', 'SingleCellExperiment', 'SummarizedExperiment','ComplexHeatmap','EnhancedVolcano'))" \
    ## clean up
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/ \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN Rscript -e "install.packages('Seurat',repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('devtools',repos='http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_github('PAIN-initiative/H5MANIPULATOR');remotes::install_github('satijalab/seurat-wrappers')"
RUN Rscript -e "remotes::install_github('PAIN-initiative/qcreporter')"
RUN git clone https://github.com/PAIN-initiative/nextflow-scRNAseq.git
RUN wget -qO- https://get.nextflow.io | bash
