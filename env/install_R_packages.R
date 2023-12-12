# First, install devtools (for installing GitHub packages) if it isn’t already installed:
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
install.packages(c('BiocManager', 'devtools','FactoMineR', 'ggforce','systemfonts', 'assertthat', 'cowplot', 'data.table', 'dplyr', 'ids', 'ggplot2', 'jsonlite', 'Matrix', 'optparse', 'GetoptLong', 'purrr','R.utils','rmarkdown','stringr','gt','plotly','tidyr','future','future.apply','rio','purrr','scCustomize','egg','DT','SoupX','reticulate','glmpca','remotes','viridis','gridExtra','plyr','circlize','naniar','XML','RCurl'),repos='http://cran.us.r-project.org')
BiocManager::install(c('rhdf5', 'SingleCellExperiment', 'SummarizedExperiment','ComplexHeatmap','EnhancedVolcano','scuttle'))
BiocManager::install("MAST")

install.packages('ggpubr',repos='http://cran.us.r-project.org')
install.packages('spam',repos='http://cran.us.r-project.org')
install.packages("harmony",repos='http://cran.us.r-project.org')
remotes::install_github("Simon-Leonard/FlexDotPlot")
install.packages('devtools',repos='http://cran.us.r-project.org')
remotes::install_github(repo = "samuel-marsh/scCustomize", ref = "release/2.0.0")
install.packages("hdf5r",repos='http://cran.us.r-project.org')

remotes::install_github('bwh-bioinformatics-hub/H5MANIPULATOR')
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github("PAIN-initiative/qcreporter")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
