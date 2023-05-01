#!/usr/bin/env Rscript

library(optparse)

options(stringsAsFactors=FALSE)


##################
# OPTION PARSING
##################

option_list <- list(
  make_option(opt_str = c("-i","--in_h5"),
              type = "character",
              default = NULL,
              help = "Input filtered_feature_bc_matrix.h5 file",
              metavar = "character"),
  make_option(opt_str = c("-l","--in_mol"),
              type = "character",
              default = NULL,
              help = "Input molecule_info.h5 file",
              metavar = "character"),
  make_option(opt_str = c("-k", "--in_key"),
              type = "character",
              default = NULL,
              help = "Input SampleSheet.csv file",
              metavar = "character"),
  make_option(opt_str = c("-j","--in_library"),
              type = "character",
              default = NULL,
              help = "Sample",
              metavar = "character")
  )

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

# Load Libraries
quiet_library <- function(...) {
  suppressPackageStartupMessages(library(...))
}
quiet_library(rhdf5)
quiet_library(H5MANIPULATOR)
quiet_library(Matrix)
quiet_library(ggplot2)
quiet_library(cowplot)
quiet_library(jsonlite)
quiet_library(tidyverse)
quiet_library(googlesheets4)
quiet_library(rio)
quiet_library(Seurat)

# Arg Parse

if(is.null(args$in_h5)) {
  in_h5 <- system.file("testdata/well1.h5", package = "H5MANIPULATOR")
  in_mol <- system.file("testdata/sample1_molecule_info.h5", package = "H5MANIPULATOR")
  in_sum <- system.file("testdata/sample1_metrics_summary.csv", package = "H5MANIPULATOR")
  in_key <- system.file("reference/SampleSheet_fallback.csv", package = "H5MANIPULATOR")
  in_sample <- "B000-P0C0W0"
} else {
  in_h5 <- args$in_h5
  in_mol <- args$in_mol
  in_key <- args$in_key
  in_library <- args$in_library
}

if(!file.exists(in_h5)) {
  stm(paste0("ERROR: Cannot find IN H5 file:", in_h5))
  stop()
}
if(!file.exists(in_mol)) {
  stm(paste0("ERROR: Cannot find IN Mol Info H5 file:", in_mol))
  stop()
}
if(!file.exists(in_key)) {
  stm(paste0("ERROR: Cannot find IN SampleSheet file:", in_key))
  stop()
}


### Load inputs

#### Load scRNA-seq Dataset
stm(paste0("Loading HDF5 from ", in_h5))
h5_list <- h5dump(in_h5)

#### Load SampleSheet
if (length(grep("https",in_key)) > 0) {
    ss <- read_sheet(in_key)
    } else if (length(grep(".xlsx",in_key)) > 0 ){
        ss <- import_list(in_key)
        ss <- ss$MetaTable_expectedCell
    } else {
        ss <- read.csv(in_key)
}
if ("Final list" %in% colnames(ss)) {
    ss <- ss %>% filter(ss$"Final list" == 1)
}

# Create Sample/Library Specific DataFrame
if (exists("in_sample")) {
    Sample_ss  <- ss %>% filter(Sample == in_sample)
} else {
    Library_ss <- ss %>% filter(Library == in_library)
}

# pull Sample/Library ID
if (exists("Sample_ss")) {
    library_id <- Sample_ss$Library
    n_cells <- length(h5_list$matrix$barcodes)
    h5_list <- set_list_path(h5_list,
                           "/matrix/observations/in_sample",
                           rep(in_sample, n_cells))
} else {
    sample_id <- Library_ss$Sample
    n_cells <- length(h5_list$matrix$barcodes)
    h5_list <- set_list_path(h5_list,
                           "/matrix/observations/in_sample",
                           rep(sample_id, n_cells))
    h5_list <- set_list_path(h5_list,
                           "/matrix/observations/in_library",
                            rep(in_library, n_cells))
}

# Pull Treatment 
if (exists("Sample_ss")) {
    Treatment_id <- Sample_ss$Treatment
    h5_list <- set_list_path(h5_list,
                            "/matrix/observations/treatment",
                            rep(Treatment_id, n_cells))

} else {
    Treatment_id <- Library_ss$Treatment
    h5_list <- set_list_path(h5_list,
                            "/matrix/observations/treatment",
                            rep(Treatment_id, n_cells))

}

# Pull Tissue 
if (exists("Sample_ss")) {
    if (length(Sample_ss$Tissue) > 0){
        Tissue <- Sample_ss$Tissue
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/tissue",
                                rep(Tissue, n_cells))
}
} else {
    if (length(Library_ss$Tissue) > 0){
        Tissue <- Library_ss$Tissue
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/tissue",
                                rep(Tissue, n_cells))
}
}

# Pull Sex 
if (exists("Sample_ss")) {
    if (length(Sample_ss$Sex) > 0){
        Sex <- Sample_ss$Sex
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/Sex",
                                rep(Sex, n_cells))
}
} else {
    if (length(Library_ss$Sex) > 0){
        Sex <- Library_ss$Sex
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/Sex",
                                rep(Sex, n_cells))
}
}

# Pull Sequencing Batch 
if (exists("Sample_ss")) {
    if (length(Sample_ss$'Sequencing Batch') > 0){
        Sequencing_Batch_id <- Sample_ss$'Sequencing Batch'
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/Sequencing_Batch_id",
                                rep(Sequencing_Batch_id, n_cells))
}
} else {
    if (length(Library_ss$'Sequencing_Batch_id') > 0){
        Sequencing_Batch_id <- Library_ss$'Sequencing_Batch_id'
        h5_list <- set_list_path(h5_list,
                                "/matrix/observations/Sequencing_Batch_id",
                                rep(Sequencing_Batch_id, n_cells))
}
}

out_h5 <- file.path(getwd(), paste0(in_library, ".h5"))
stm(paste0("OUT H5 file         : ", out_h5))

#### Read molecule info to get read counts per cell
stm(paste0("Assembling Read Counts per Cell from ", in_mol))
bc <- sub("-1","",h5_list$matrix$barcodes)
bc_counts <- data.table(mol_idx = h5read(in_mol, "/barcode_idx"),
                        umi_count = h5read(in_mol, "/count"))
bc_sums <- bc_counts[, .(n_reads = sum(umi_count)), by = mol_idx]
rm(bc_counts)
mol_bc <- h5read(in_mol, "/barcodes")
bc_sums$cell_barcode <- mol_bc[bc_sums$mol_idx + 1]
rm(mol_bc)
bc_sums <- bc_sums[,.(cell_barcode, n_reads)]
n_reads <- bc_sums$n_reads[match(bc, bc_sums$cell_barcode)]
n_reads[is.na(n_reads)] <- 0
h5_list <- set_list_path(h5_list,
                           "/matrix/observations/n_reads",
                           n_reads)

### Assemble data

#### Split matrices if more than one feature type is present

h5_list <- h5_list_convert_to_dgCMatrix(h5_list, target = "matrix")
feature_types <- unique(h5_list$matrix$features$feature_type)
if(sum(feature_types != "Gene Expression") > 0) {
  stm("Separating non-Gene Expression data to additional matrices")
  
  mat <- h5_list$matrix_dgCMatrix
  
  feature_df <- as.data.frame(h5_list$matrix$features)
  
  h5_list$matrix_dgCMatrix <- mat[feature_df$feature_type == "Gene Expression",]
  h5_list$matrix$features <- as.list(feature_df[feature_df$feature_type == "Gene Expression",])
}

#### Compute N UMIs and N Genes per cell
stm("Computing UMI and Gene Counts per Cell")
h5_list <- set_list_path(h5_list,
                         "/matrix/observations/n_umis",
                         unname(colSums(h5_list$matrix_dgCMatrix)))
h5_list <- set_list_path(h5_list,
                         "/matrix/observations/n_genes",
                         unname(colSums(h5_list$matrix_dgCMatrix > 0)))
h5_list <- h5_list_convert_from_dgCMatrix(h5_list, target = "matrix")

#### Add cell ids
stm("Adding Cell UUIDs and Names")
h5_list <- add_cell_ids(h5_list,
                        add_uuid = TRUE,
                        replace_barcode = TRUE,
                        retain_original_barcode = TRUE,
                        add_name = TRUE)

#### Add chrM gene counts
stm("Adding chrM count metadata")
if ((unique(h5_list$matrix$features$genome) == "GRCh38") == TRUE){
    h5_list <- h5_list_add_mito_umis(h5_list)
    } else {
        so <- read_h5_seurat(in_h5)
        all_genes <- h5_list$matrix$features$name
        mito_genes <- grep("^MT:|MT-", all_genes,ignore.case=TRUE,value=TRUE)
        total_counts_per_cell <- colSums(so@assays$RNA@counts)
        percent_mito <- colSums(so@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
        h5_list <- H5MANIPULATOR::set_list_path(h5_list,
                                     "/matrix/observations/n_mito_umis",
                                     percent_mito)
}


#### Add Sample Metrics
if (exists("in_sum")) {
    sample_metrics <- read_tenx_metrics(in_sum)
    sample_metrics <- as.list(sample_metrics)
    h5_list <- set_list_path(h5_list,
                            "/sample",
                            sample_metrics)
    h5_list <- set_list_path(h5_list,
                            "/sample/sample_id",
                            in_sample)
} else {
    estimated_number_of_cells   <- Library_ss$"Estimated Number of Cells"
    h5_list <- set_list_path(h5_list,
                            "/sample/estimated_number_of_cells",
                            estimated_number_of_cells)
    fraction_reads_in_cells                         <- Library_ss$"Fraction Reads in Cells"
    h5_list <- set_list_path(h5_list,
                            "/sample/fraction_reads_in_cells",
                            fraction_reads_in_cells)
    mean_reads_per_cell                             <- Library_ss$"Mean Reads per Cell"
    h5_list <- set_list_path(h5_list,
                            "/sample/mean_reads_per_cell",
                            mean_reads_per_cell)
    median_genes_per_cell                           <- Library_ss$"Median Genes per Cell"
    h5_list <- set_list_path(h5_list,
                            "/sample/median_genes_per_cell",
                            median_genes_per_cell)
    median_umi_counts_per_cell                      <- Library_ss$"Estimated Number of Cells"
    h5_list <- set_list_path(h5_list,
                            "/sample/median_umi_counts_per_cell",
                            median_genes_per_cell)
    number_of_reads                                 <- Library_ss$"Number of Reads"
    h5_list <- set_list_path(h5_list,
                            "/sample/number_of_reads",
                            number_of_reads)
    q30_bases_in_barcode                            <- Library_ss$"Q30 Bases in Barcode"
    h5_list <- set_list_path(h5_list,
                            "/sample/q30_bases_in_barcode",
                            q30_bases_in_barcode)
    q30_bases_in_rna_read                           <- Library_ss$"Q30 Bases in RNA Read"
    h5_list <- set_list_path(h5_list,
                            "/sample/q30_bases_in_rna_read",
                            q30_bases_in_rna_read)
    q30_bases_in_umi                                <- Library_ss$"Q30 Bases in UMI"
    h5_list <- set_list_path(h5_list,
                            "/sample/q30_bases_in_umi",
                            q30_bases_in_umi)
    reads_mapped_antisense_to_gene                  <- Library_ss$"Reads Mapped Antisense to Gene"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_antisense_to_gene",
                            reads_mapped_antisense_to_gene)
    reads_mapped_confidently_to_exonic_regions      <- Library_ss$"Reads Mapped Confidently to Exonic Regions"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_confidently_to_exonic_regions",
                            reads_mapped_confidently_to_exonic_regions)
    reads_mapped_confidently_to_genome              <- Library_ss$"Reads Mapped Confidently to Genome"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_confidently_to_genome",
                            reads_mapped_confidently_to_genome)
    reads_mapped_confidently_to_intergenic_regions  <- Library_ss$"Reads Mapped Confidently to Intergenic Regions"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_confidently_to_intergenic_regions",
                            reads_mapped_confidently_to_intergenic_regions)
    reads_mapped_confidently_to_intronic_regions    <- Library_ss$"Reads Mapped Confidently to Intronic Regions"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_confidently_to_intronic_regions",
                            reads_mapped_confidently_to_intronic_regions)
    reads_mapped_confidently_to_transcriptome       <- Library_ss$"Reads Mapped Confidently to Transcriptome"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_confidently_to_transcriptome",
                            reads_mapped_confidently_to_transcriptome)
    reads_mapped_to_genome                          <- Library_ss$"Reads Mapped Confidently to Genome"
    h5_list <- set_list_path(h5_list,
                            "/sample/reads_mapped_to_genome",
                            reads_mapped_to_genome)
    sample_id                                       <- Library_ss$"Sample"
    h5_list <- set_list_path(h5_list,
                         "/sample/sample_id",
                         sample_id)
    sequencing_saturation                           <- Library_ss$"Sequencing Saturation"
    h5_list <- set_list_path(h5_list,
                         "/sample/sequencing_saturation",
                         sequencing_saturation)
    total_genes_detected                            <- Library_ss$"Total Genes Detected"
    h5_list <- set_list_path(h5_list,
                         "/sample/total_genes_detected",
                         total_genes_detected)
    valid_barcodes                                  <- Library_ss$"Valid Barcodes"
    h5_list <- set_list_path(h5_list,
                         "/sample/valid_barcodes",
                         valid_barcodes)
}

#### Write HDF5 files
stm(paste0("Writing HDF5 to ", out_h5))
write_h5_list(h5_list,
              h5_file = out_h5,
              overwrite = TRUE)
h5closeAll()

