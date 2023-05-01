FROM ubuntu

## Install apt dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    apt-utils \
    gfortran \
    liblapack-dev \
    liblapack3 \
    libopenblas-base \
    libopenblas-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    pandoc \
    libhdf5-dev \
    git

## clean up
RUN apt-get clean \
  && rm -rf /var/lib/apt/lists/ \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 8B57C5C2836F4BEB
# RUN apt-get update
## Install R packages
RUN export R_HOME=/usr/lib/R
RUN R -e 'install.packages(c("devtools"))'

## CRAN dependencies
RUN R -e 'install.packages(c("jsonlite","rmarkdown","optparse","viridis"))'

RUN R -e 'install.packages(c("devtools","BiocManager"));BiocManager::install()'

## BioConductor dependencies
RUN R -e 'BiocManager::install(c("rhdf5","GenomicRanges"))'



## Install Required GITHUB packages
COPY auth_token /tmp/auth_token
RUN export GITHUB_PAT=$(cat /tmp/auth_token) \
   && R -e    'auth_token = Sys.getenv("GITHUB_PAT")); devtools::install_github("bwh-bioinformatics-hub/H5MANIPULATOR", auth_token =   Sys.getenv("GITHUB_PAT")); devtools::install_github("acicalo2/qcreporter", auth_token = Sys.getenv("GITHUB_PAT"))' \
  && git clone  https://aifi-gitops:$GITHUB_PAT@github.com/bwh-bioinformatics-hub/qcreporter.git \
  && rm -rf /tmp/downloaded_packages /tmp/*.rds /tmp/auth_token 

## Pipeline package requirements
RUN R -e 'install.packages(c("rmarkdown","optparse"))'


