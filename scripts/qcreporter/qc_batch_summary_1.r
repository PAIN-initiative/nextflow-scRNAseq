#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(opt_str = c("-e","--experiment_id"),
              type = "character",
              default = NULL,
              help = "Sample identifier",
              metavar = "character"),
  make_option(opt_str = c("-m","--in_method"),
              type = "character",
              default = NULL,
              help = "Input batch pipeline modality string",
              metavar = "character"),
  make_option(opt_str = c("-i","--in_dir"),
              type = "character",
              default = NULL,
              help = "Input directory containing h5 and json files",
              metavar = "character"),
  make_option(opt_str = c("-z","--cellrangers_dir"),
              type = "character",
              default = NULL,
              help = "cellrangers out directory",
              metavar = "character"),
  make_option(opt_str = c("-f","--refdir"),
              type = "character",
              default = NULL,
              help = "reference directory",
              metavar = "character"),
  make_option(opt_str = c("-k","--in_key"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
              metavar = "character"),
  make_option(opt_str = c("-d","--out_dir"),
              type = "character",
              default = NULL,
              help = "Output Directory",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_html"),
              type = "character",
              default = NULL,
              help = "Output HTML run summary file",
              metavar = "character"),
  make_option(opt_str = c("-j","--resolution"),
              type = "character",
              default = NULL,
              help = "Resolution choice to set identity of Seurat Object",
              metavar = "character"),
  make_option(opt_str = c("-c","--percent_mito"),
              type = "integer",
              default = 10,
              help = "Percent Mito Filter Value",
              metavar = "integer"),
  make_option(opt_str = c("-a","--percent_ribo"),
              type = "integer",
              default = 20,
              help = "Percent Ribo Filter Value",
              metavar = "integer"),
  make_option(opt_str = c("-b","--filter_MALAT"),
              type = "character",
              default = FALSE,
              help = "Filter MALAT Gene TRUE/FALSE",
              metavar = "character"),
  make_option(opt_str = c("-u","--filter_MITO"),
              type = "character",
              default = FALSE,
              help = "Filter Mito Genes TRUE/FALSE",
              metavar = "character"),
  make_option(opt_str = c("-l","--species"),
              type = "character",
              default = NULL,
              help = "Species of Experiment",
              metavar = "character"),
  make_option(opt_str = c("-q","--filter_RIBO"),
              type = "character",
              default = FALSE,
              help = "Filter Ribo Genes TRUE/FALSE",
              metavar = "character"))

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)
# Load Libraries
quiet_library <- function(...) {
  suppressPackageStartupMessages(library(...))
}
quiet_library(optparse)        # dependency of H5MANIPULATOR
quiet_library(qcreporter)        # dependency of H5MANIPULATOR
quiet_library(Matrix)        # dependency of H5MANIPULATOR
quiet_library(rhdf5)         # dependency of H5MANIPULATOR
quiet_library(H5MANIPULATOR)    
quiet_library(ggplot2)
quiet_library(stringr)       
quiet_library(dplyr)         # data wrangling
quiet_library(cowplot)       # arranging multiple plots
quiet_library(gt)            # formatted table output
quiet_library(plotly)        # interactive plots
quiet_library(tidyr)         # data wrangling
quiet_library(Seurat)        # batch umap creation
quiet_library(future)        # multi-threading for batch umap creation
quiet_library(future.apply)  # multi-threading for batch umap creation
quiet_library(rio)
quiet_library(purrr)
quiet_library(egg)
quiet_library(DT)
quiet_library(SoupX)
quiet_library(reticulate)
quiet_library(glmpca)
quiet_library(SeuratWrappers)
quiet_library(viridis)
quiet_library(qs)
quiet_library(gridExtra)
quiet_library(plyr)
quiet_library(ggpubr)
quiet_library(XML)
quiet_library(RCurl)
quiet_library(DoubletFinder)
quiet_library(scuttle)


# Arg Parse

experiment_id   <- args$experiment_id
in_method 		  <- args$in_method
in_dir			    <- args$in_dir
cellrangers_dir <- args$cellrangers_dir
refdir 			    <- args$refdir
in_key 			    <- args$in_key
out_dir 		    <- args$out_dir
out_html 		    <- args$out_html
resolution 		  <- args$resolution
percent_mito 	  <- args$percent_mito
percent_ribo 	  <- args$percent_ribo
filter_MALAT 	  <- args$filter_MALAT
filter_MITO 	  <- args$filter_MITO
species 		    <- args$species
filter_RIBO 	  <- args$filter_RIBO
# Load Samplesheet
stm("Reading in Marker Genes")

if (length(grep("https",in_key)) > 0) {
    ss <- read_sheet(in_key)
    } else if (length(grep(".xlsx",in_key)) > 0 ){
        ss <- import_list(in_key)
    } else {
        ss <- read.csv(in_key)
}

metatable <- ss$MetaTable_expectedCell
if ("Final list" %in% colnames(metatable)) {
    metatable <- metatable %>% filter(metatable$"Final list" == 1)
}
samples <- unique(metatable$Sample)
treatments <- metatable$Treatment      
              
# Format Marker Genes
MarkerGenes <- ss$MarkerGenes
colnames(MarkerGenes) <- c('marker_gene','cell_type','ref')
markers <- MarkerGenes$'marker_gene'

# read in gene name table
geneTable <- read.csv(paste0(refdir, "geneAnnotationTable.csv"), header = T, row.names = 1)

# scCustermize code updated by developer, https://github.com/samuel-marsh/scCustomize/issues/27
# Clustered_DotPlot_relabel function from developer

Clustered_DotPlot_relabel <- function(
  seurat_object,
  features,
  new_row_labels,
  colors_use_exp = viridis_plasma_dark_high,
  exp_color_min = -2,
  exp_color_middle = NULL,
  exp_color_max = 2,
  print_exp_quantiles = FALSE,
  colors_use_idents = NULL,
  x_lab_rotate = TRUE,
  k = 1,
  row_km_repeats = 1000,
  column_km_repeats = 1000,
  row_label_size = 8,
  raster = FALSE,
  plot_km_elbow = TRUE,
  elbow_kmax = NULL,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  show_parent_dend_line = TRUE,
  ggplot_default_colors = FALSE,
  seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- PackageCheck("ComplexHeatmap", error = FALSE)
  if (!ComplexHeatmap_check[1]) {
    stop(
      "Please install the ComplexHeatmap package to use Clustered_DotPlot",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('ComplexHeatmap')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  
  # Check Seurat
  scCustomize:::Is_Seurat(seurat_object = seurat_object)
  
  # Check unique features
  features_unique <- unique(x = features)
  
  if (length(x = features_unique) != length(x = features)) {
    warning("Feature list contains duplicates, making unique.")
  }
  
  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    stop("The value for 'exp_color_min': ", exp_color_min, ", must be less than the value for 'exp_color_max': ", exp_color_max, ".")
  }
  
  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = features_unique, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)
  
  data <- seurat_plot$data
  
  # Get expression data
  exp_mat <- data %>%
    dplyr::select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  
  row.names(x = exp_mat) <- exp_mat$features.plot
  
  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    warning("The following features were removed as there is no scaled expression present in subset (`idents`) of object provided: ", glue_collapse_scCustom(input_string = excluded_features, and = TRUE), ".")
    
    # Extract good features
    good_features <- rownames(exp_mat)
    
    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(features.plot %in% good_features)
  }
  
  exp_mat <- exp_mat[,-1] %>%
    as.matrix()
  
  # Get percent expressed data
  percent_mat <- data %>%
    dplyr::select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()
  
  row.names(x = percent_mat) <- percent_mat$features.plot
  
  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(features.plot %in% good_features)
  }
  
  percent_mat <- percent_mat[,-1] %>%
    as.matrix()
  
  # print quantiles
  if (print_exp_quantiles) {
    message("Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }
  
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    if (is.null(x = colors_use_idents)) {
      colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }
  
  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(exp_mat)
  
  identity_colors <- DiscretePalette_scCustomize(num_colors = length(Identity), palette = "polychrome", shuffle_pal = F)
  names(identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)
  
  # Create identity annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                 col =  identity_colors_list,
                                                 na_col = "grey",
                                                 name = "Identity"
  )
  
  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- scCustomize:::Middle_Number(min = exp_color_min, max = exp_color_max)
  }
  
  palette_length <- length(colors_use_exp)
  palette_middle <- scCustomize:::Middle_Number(min = 0, max = palette_length)
  
  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
  
  # Calculate and plot Elbow
  if (plot_km_elbow) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      warning("The value provided for 'elbow_kmax' is too large.  Changing to (length(x = features)-1): ", elbow_kmax)
    }
    
    # if elbow_kmax is NULL set value based on input feature list
    if (is.null(x = elbow_kmax)) {
      # set to (length(x = features)-1) if less than 21 features OR to 20 if greater than 21 features
      if (nrow(x = exp_mat) > 21) {
        elbow_kmax <- 20
      } else {
        elbow_kmax <- nrow(x = exp_mat) - 1
      }
    }
    
    km_elbow_plot <- scCustomize:::kMeans_Elbow(data = exp_mat, k_max = elbow_kmax)
  }
  
  # prep heatmap
  if (raster) {
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                  gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
    }
  } else {
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                  gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
    }
  }
  
  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(labels = c(0.25,0.5,0.75,1), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black")))
    )
  )
  
  # Set x label roration
  if (is.numeric(x = x_lab_rotate)) {
    x_lab_rotate <- x_lab_rotate
  } else if (isTRUE(x = x_lab_rotate)) {
    x_lab_rotate <- 45
  } else {
    x_lab_rotate <- 0
  }
  
  # Create Plot
  set.seed(seed = seed)
  if (raster) {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                layer_fun = layer_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  } else {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                cell_fun = cell_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  }
  
  # Add pt.size legend & return plots
  if (plot_km_elbow) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot + rowAnnotation(rn= anno_text(new_row_labels)), annotation_legend_list = lgd_list))
}

# Step 2: Pre-processing
# Remove ambient RNA by SoupX
data.10x = list()
for (sample in samples){
  filt.matrix <- Read10X_h5(paste0(cellrangers_dir, sample, "/outs/filtered_feature_bc_matrix.h5"), use.names = T)
  raw.matrix <- Read10X_h5(paste0(cellrangers_dir, sample, "/outs/raw_feature_bc_matrix.h5"), use.names = T)
  srat <- CreateSeuratObject(counts = filt.matrix)
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)
  srat <- SCTransform(srat, verbose = F)
  srat <- RunPCA(srat, verbose = F)
  srat <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat <- FindClusters(srat, verbose = T)
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel <- autoEstCont(soup.channel)
  data.10x[[sample]] <- adjustCounts(soup.channel, roundToInt = T)
}

# Create Seurat object after SoupX
scrna.list = list()
for (sample in samples) {
    scrna.list[[sample]] = CreateSeuratObject(counts = data.10x[[sample]], min.cells=3, project=sample)
}

# Remove raw data to save memory
rm(data.10x)
                        
# add treatment
for(i in 1:length(samples)){
  sample=samples[i]; treatment=treatments[i];
  scrna.list[[sample]]$treatment <- treatment
}

# add sample name
for(i in 1:length(samples)){
  sample=samples[i]; sample_id=samples[i];
  scrna.list[[sample]]$sample_id <- sample_id
}

# Add percent.mt and percent.rb to cell level metadata
for (sample in samples) {
  scrna.list[[sample]][["percent.mito"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^MT:|MT-|mt:|mt-") 
  scrna.list[[sample]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^RP[LS]|Rp[LS]")
}
# merge list of prefiltered Seurat
scrna.combined_prefilter <- Merge_Seurat_List(scrna.list,
                                              add.cell.ids = NULL,
                                              merge.data = TRUE,
                                              project = "sample_id")
# metadata variable
metadata_prefilter <- scrna.combined_prefilter@meta.data
# Save Prefilter
saveRDS(scrna.combined_prefilter, paste0(out_dir, "scrna.combined_prefilter.seurat.", projectName, ".rds"))
                   
# Filtered Cells with 3SD of mean nCount and nFeature, percent of mito
  
# Filtered cells with 3SD of mean nCount and nFeature, percent of mito
qc_cutoff = 3
mito_cutoff = percent_mito
ribo_cutoff = percent_ribo
for (sample in samples){
  mean.nCount <- mean(scrna.list[[sample]]@meta.data$nCount_RNA)
  sd.nCount <- sd(scrna.list[[sample]]@meta.data$nCount_RNA)
  mean.nFeature <- mean(scrna.list[[sample]]@meta.data$nFeature_RNA)
  sd.nFeature <- sd(scrna.list[[sample]]@meta.data$nFeature_RNA)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset = nCount_RNA > mean.nCount - qc_cutoff*sd.nCount & nCount_RNA < mean.nCount + qc_cutoff*sd.nCount & nFeature_RNA > mean.nFeature - qc_cutoff*sd.nFeature & nFeature_RNA < mean.nFeature + qc_cutoff*sd.nFeature & percent.mito < mito_cutoff & percent.ribo < ribo_cutoff)
}
  
# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
C = list()
most_expressed = list()
par(mar = c(4, 8, 2, 1))
for (sample in samples){
C[[sample]] <- scrna.list[[sample]]@assays$RNA@counts
C[[sample]] <- Matrix::t(Matrix::t(C[[sample]])/Matrix::colSums(C[[sample]])) * 100
most_expressed[[sample]] <- order(apply(C[[sample]], 1, median), decreasing = T)[20:1]
}
most_expressed_plots = list()
# most_expressed_plots[[sample]] <- 
for (sample in samples){
    pdf(paste0(out_dir,sample,"_most_expressed_genes.pdf"))
    boxplot(as.matrix(t(C[[sample]][most_expressed[[sample]], ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
    dev.off()
}

    # Filter MALAT1
if (filter_MALAT == TRUE){
  for (sample in samples){
  scrna.list[[sample]] <- scrna.list[[sample]][!grepl("MALAT1", rownames(scrna.list[[sample]])), ]
  }
}
# Filter Mitocondrial
if (filter_MITO == TRUE){
  for (sample in samples){
  scrna.list[[sample]] <- scrna.list[[sample]][!grepl("^MT:|MT-|mt:|mt-", rownames(scrna.list[[sample]])), ]
  }
}
# Filter Ribosomal gene (optional if that is a problem on your data) data.filt
if (filter_RIBO == TRUE){
  for (sample in samples){
  scrna.list[[sample]] <- scrna.list[[sample]][!grepl("^RP[LS]|Rp[LS]", rownames(scrna.list[[sample]])), ]
  }
}
### Doublet Finder
# grab number of cells
n_cells = list()
for (sample in samples){
n_cells[[sample]] <- length(colnames(scrna.list[[sample]]))
}
multiplet_rate = list()
for (sample in samples){
    if (n_cells[[sample]] <= 500){
        multiplet_rate[[sample]] = 0.004
    }else if (n_cells[[sample]] > 500 & n_cells[[sample]] < 2000){
             multiplet_rate[[sample]] = 0.008 
    } else if (n_cells[[sample]] >= 2000 & n_cells[[sample]] < 3000){
      multiplet_rate[[sample]] = 0.016 
    } else if (n_cells[[sample]] >= 3000 & n_cells[[sample]] < 4000){
      multiplet_rate[[sample]] = 0.024
    } else if (n_cells[[sample]] >= 4000 & n_cells[[sample]] < 5000){
      multiplet_rate[[sample]] = 0.032
    } else if (n_cells[[sample]] >= 5000 & n_cells[[sample]] < 6000){
      multiplet_rate[[sample]] = 0.040
    } else if (n_cells[[sample]] >= 6000 & n_cells[[sample]] < 7000){
      multiplet_rate[[sample]] = 0.048
    } else if (n_cells[[sample]] >= 7000 & n_cells[[sample]] < 8000){
      multiplet_rate[[sample]] = 0.056
    } else if (n_cells[[sample]] >= 8000 & n_cells[[sample]] < 9000){
      multiplet_rate[[sample]] = 0.064
    } else if (n_cells[[sample]] >= 9000 & n_cells[[sample]] < 10000){
      multiplet_rate[[sample]] = 0.072
    }else{
     multiplet_rate[[sample]] = 0.080
    }
}
nExp = list()
for (sample in samples){
scrna.list[[sample]] <- scrna.list[[sample]] %>% NormalizeData()
scrna.list[[sample]] = FindVariableFeatures(scrna.list[[sample]], verbose = F)
scrna.list[[sample]] = ScaleData(scrna.list[[sample]],verbose = F)
scrna.list[[sample]] = RunPCA(scrna.list[[sample]], verbose = F, npcs = 20)
scrna.list[[sample]] = RunUMAP(scrna.list[[sample]], dims = 1:10, verbose = F)
nExp[[sample]] <- round(ncol(scrna.list[[sample]]) * multiplet_rate[[sample]])  # expected doublets
scrna.list[[sample]] <- suppressMessages(doubletFinder_v3(scrna.list[[sample]], pN = 0.25, pK = 0.09, nExp = nExp[[sample]], PCs = 1:10))
}


DF.name = list()
for (sample in samples){
# name of the DF prediction can change, so extract the correct column name.
DF.name[[sample]] = colnames(scrna.list[[sample]]@meta.data)[grepl("^DF.classification", colnames(scrna.list[[sample]]@meta.data))]
}
# Plot the Doublet Finder results
UMAP_plots <- list()
for (sample in samples){
    UMAP_plots[[sample]] <-  cowplot::plot_grid(ncol = 2, DimPlot(scrna.list[[sample]], group.by = "sample_id") + NoAxes(),
    DimPlot(scrna.list[[sample]], group.by = DF.name[[sample]]) + NoAxes())
}

for (sample in samples){
    pdf(paste0(out_dir,sample,"_DoubletFinder_UMAP_Plot.pdf"),width=15,height=15,onefile=TRUE)
    print(UMAP_plots[[sample]])
    dev.off()
}

VlnPlots = list()
for (sample in samples){
VlnPlots[[sample]] <- VlnPlot(scrna.list[[sample]], features = "nFeature_RNA", group.by = DF.name[[sample]], pt.size = 0.1)
}

for (sample in samples){
    pdf(paste0(out_dir,sample,"_DoubletFinder_VlnPlot.pdf"),width=15,height=15,onefile=TRUE)
    print(VlnPlots[[sample]])
    dev.off()
}

# Remove the Doublet Cells
cells.use = list()
for (sample in samples){
     cells.use[[sample]] <- colnames(scrna.list[[sample]])[which(scrna.list[[sample]][[]][DF.name[[sample]]] == "Singlet")]
     scrna.list[[sample]] <- subset(scrna.list[[sample]], cells = cells.use[[sample]])
}

so_postfilter <- Merge_Seurat_List(scrna.list,
                                              add.cell.ids = NULL,
                                              merge.data = TRUE,
                                              project = "sample_id")
metadata_postfilter <- so_postfilter@meta.data

# mito
prefilter_qc_mito <- metadata_prefilter %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mito)) + 
geom_point() + 
scale_colour_gradient(low = "gray90", high = "black") +
stat_smooth(method=lm) +
scale_x_log10() +
scale_y_log10() + 
theme_classic() +
geom_vline(xintercept = 500) +
geom_hline(yintercept = 300) +
facet_wrap(~sample_id) + 
ggtitle("Prefilter Mitochondrial Percentage")

postfilter_qc_mito <- metadata_postfilter %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mito)) + 
geom_point() + 
scale_colour_gradient(low = "gray90", high = "black") +
stat_smooth(method=lm) +
scale_x_log10() +
scale_y_log10() + 
theme_classic() +
geom_vline(xintercept = 500) +
geom_hline(yintercept = 300) +
facet_wrap(~sample_id) + 
ggtitle("Postfilter Mitochondrial Percentage")

pdf(paste0(out_dir,"Mitochondrial_Percentage.pdf"))
ggarrange(prefilter_qc_mito,postfilter_qc_mito,ncol=1,nrow=2,widths = 25,heights = 30)
dev.off()

# ribo
prefilter_qc_ribo <- metadata_prefilter %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.ribo)) + 
geom_point() + 
scale_colour_gradient(low = "gray90", high = "black") +
stat_smooth(method=lm) +
scale_x_log10() +
scale_y_log10() + 
theme_classic() +
geom_vline(xintercept = 500) +
geom_hline(yintercept = 300) +
facet_wrap(~sample_id) +
ggtitle("Prefilter Ribosomal Percentage")

postfilter_qc_ribo <- metadata_postfilter %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.ribo)) + 
geom_point() + 
scale_colour_gradient(low = "gray90", high = "black") +
stat_smooth(method=lm) +
scale_x_log10() +
scale_y_log10() + 
theme_classic() +
geom_vline(xintercept = 500) +
geom_hline(yintercept = 300) +
facet_wrap(~sample_id) +
ggtitle("Postfilter Ribosomal Percentage")

pdf(paste0(out_dir,"Ribosomal_Percentage.pdf"))
ggarrange(prefilter_qc_ribo,postfilter_qc_ribo,ncol=1,nrow=2,widths = 10,heights = 8)
dev.off()

UMIs_transcripts_per_cell_prefilter <- metadata_prefilter %>% 
  	ggplot(aes(color=sample_id, x=nCount_RNA, fill= sample_id)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500) + 
    ggtitle("Prefilter UMIs Transcripts Per Cell")

# UMIs_transcripts_per_cell
UMIs_transcripts_per_cell_postfilter <- metadata_postfilter %>% 
  	ggplot(aes(color=sample_id, x=nCount_RNA, fill= sample_id)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500) + 
    ggtitle("Postfilter UMIs Transcripts Per Cell")

pdf(paste0(out_dir,"UMIs_transcripts_per_cell.pdf"))
ggarrange(UMIs_transcripts_per_cell_prefilter,UMIs_transcripts_per_cell_postfilter,ncol=1,nrow=2,widths = 10,heights = 8)
dev.off()

# Genes Per Cell
genes_per_cell_pre <- metadata_prefilter %>% 
  	ggplot(aes(color=sample_id, x=nCount_RNA, fill=sample_id )) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300) + 
    ggtitle("Prefilter Genes Per Cell")

genes_per_cell_post <- metadata_postfilter %>% 
  	ggplot(aes(color=sample_id, x=nCount_RNA, fill=sample_id )) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300) + 
    ggtitle("Postfilter Genes Per Cell")

pdf(paste0(out_dir,"Genes_Per_Cell.pdf"))
ggarrange(genes_per_cell_pre,genes_per_cell_post,ncol=1,nrow=2,widths = 10,heights = 8)
dev.off()

# Integration
# normalize and identify variable features for each dataset independently
scrna.list <- SplitObject(so_postfilter, split.by = "sample_id")
scrna.list <- lapply(X = scrna.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scrna.list)
# Perform integration
scrna.anchors <- FindIntegrationAnchors(object.list = scrna.list, anchor.features = features)
scrna.combined <- IntegrateData(anchorset = scrna.anchors)
rm(scrna.anchors)
# Perform an integrated analysis
DefaultAssay(scrna.combined) <- "integrated"
scrna.combined <- ScaleData(scrna.combined, verbose = FALSE)
scrna.combined <- RunPCA(scrna.combined, npcs = 15, verbose = FALSE)

# Linear Dimensionality Reduction
pdf(paste0(out_dir, "ElbowPlot.pdf"))
p1 <- ElbowPlot(scrna.combined) + ggtitle("Integrated") + theme(aspect.ratio=5/10) + theme(plot.margin = unit(c(6, 6, 6, 6), "cm"))
print(p1 + coord_fixed())
dev.off()

### UMAP - Resolutions
# Continue on analysis
scrna.combined <- RunUMAP(scrna.combined, reduction = "pca", dims = 1:8,verbose=FALSE)
scrna.combined <- FindNeighbors(scrna.combined, reduction = "pca", dims = 1:8,verbose=FALSE)
scrna.combined <- FindClusters(scrna.combined, resolution = seq(0.5,3,by=0.1),verbose=FALSE)
res_0.5 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.0.5", label = TRUE)
pdf(paste0(out_dir, "res_0.5.pdf"))
print(res_0.5)
dev.off()
res_0.5
res_0.7 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.0.7", label = TRUE)
pdf(paste0(out_dir, "res_0.7.pdf"))
print(res_0.7)
dev.off()
res_1 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.1", label = TRUE)
pdf(paste0(out_dir, "res_1.pdf"))
print(res_1)
dev.off()
res_1
res_1.2 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.1.2", label = TRUE)
pdf(paste0(out_dir, "res_1.2.pdf"))
print(res_1.2)
dev.off()
res_1.4 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.1.4", label = TRUE)
pdf(paste0(out_dir, "res_1.4.pdf"))
print(res_1.4) 
dev.off()
res_2 <- DimPlot(scrna.combined, group.by = "integrated_snn_res.2", label = TRUE)
pdf(paste0(out_dir, "res_2.pdf"))
print(res_2)
dev.off()
res_2
meta <- scrna.combined@meta.data
Idents(scrna.combined) <- resolution
rm(scrna.list)

p1 <- DimPlot(scrna.combined, reduction = "umap", group.by = "sample_id")
pdf(file=paste0(out_dir, "combined.umap.colorBySample.pdf"))
print(p1 + coord_fixed())
dev.off()
print(p1 + coord_fixed())
p2 <- DimPlot(scrna.combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf(paste0(out_dir, "combined.umap.colorByCluster.pdf"))
print(p2 + coord_fixed())
dev.off()
print(p2 + coord_fixed())
p4 <- DimPlot(scrna.combined, reduction = "umap", group.by = "treatment")
pdf(file=paste0(out_dir, "combined.umap.colorByTreatment.pdf"))
print(p4 + coord_fixed())
dev.off() 
print(p4 + coord_fixed())
saveRDS(scrna.combined, paste0(out_dir, "scrna.combined_postfilter.seurat.", projectName, ".rds"))

# Find Markers
DefaultAssay(scrna.combined) <- "RNA"
all.genes <- rownames(scrna.combined)
scrna.combined <- ScaleData(scrna.combined,features = all.genes,verbose = FALSE)
#nk.markers <- FindConservedMarkers(scrna.combined, ident.1 = 0, grouping.var = "sample_id", verbose = FALSE)
# FindAllMarkers
scrna.markers <- FindAllMarkers(scrna.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,verbose=FALSE,test.use = "MAST")
names(scrna.markers)[names(scrna.markers) == "gene"] <- "geneSymbol"
scrna.markers$geneID <- geneTable$geneID[match(scrna.markers$geneSymbol, geneTable$geneSymbol)]
write.table(scrna.markers, paste0(out_dir, "FindAllMarkers.clusters.xls"), sep = "\t", row.names = F)
# top 10
top10<- scrna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, paste0(out_dir, "FindAllMarkers.clusters.top10.xls"), sep = "\t", col.names = NA)

# top 25
top25 <- scrna.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
write.table(top25, paste0(out_dir, "FindAllMarkers.clusters.top25.xls"), sep = "\t", col.names = NA)

# top 50
top50 <- scrna.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.table(top50, paste0(out_dir, "FindAllMarkers.clusters.top50.xls"), sep = "\t", col.names = NA)

# Step 6: Top 3 identified genes, feature plot, dotplot
topN <- Extract_Top_Markers(scrna.markers, num_genes = 25, named_vector = FALSE, make_unique = TRUE, gene_column = "geneSymbol")

scrna.markers_sig <- scrna.markers[scrna.markers$p_val_adj <= 0.05,]            
# Create Interactive Table
if (species == "Mus musculus"){
tibble::as_tibble(unique(scrna.markers_sig)) %>% dplyr::arrange(p_val_adj) -> scrna.markers_sig
scrna.markers_sig = scrna.markers_sig %>% dplyr::select(cluster, geneSymbol, p_val, avg_log2FC,
                                             p_val_adj,pct.1,pct.2,geneID)
final_df <- data.frame(           "cluster"    = scrna.markers_sig$cluster,
                                  "gene"       = scrna.markers_sig$geneSymbol,
                                  "p_val"      = scrna.markers_sig$p_val,
                                  "avg_log2FC" = scrna.markers_sig$avg_log2FC,
                                  "p_val_adj"  = scrna.markers_sig$p_val_adj,
                                  "pct.1"      = scrna.markers_sig$pct.1,
                                  "pct.2"      = scrna.markers_sig$pct.2,
                                  "MGI_ID" = paste0("<a  href='https://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=",scrna.markers_sig$geneID,"'>", scrna.markers_sig$geneID,"</a>"))

dtable <- DT::datatable(final_df,class = 'cell-border stripe',rownames=F,filter='top',
              editable = TRUE, extensions = 'Buttons', options = list(
                dom = 'Bfrtip',
                buttons = c('copy','csv','excel','pdf','print')
              ),escape = FALSE)
}

saveWidget(dtable,"FindMarkers.html")

if (species == "Drosophila Melanogaster"){
as_tibble(unique(scrna.markers_sig)) %>% arrange(p_val_adj) -> scrna.markers_sig
scrna.markers_sig = scrna.markers_sig %>% dplyr::select(cluster, geneSymbol, p_val, avg_log2FC,
                                             p_val_adj,pct.1,pct.2,geneID)
final_df <- data.frame(           "cluster"    = scrna.markers_sig$cluster,
                                  "gene"       = scrna.markers_sig$geneSymbol,
                                  "p_val"      = scrna.markers_sig$p_val,
                                  "avg_log2FC" = scrna.markers_sig$avg_log2FC,
                                  "p_val_adj"  = scrna.markers_sig$p_val_adj,
                                  "pct.1"      = scrna.markers_sig$pct.1,
                                  "pct.2"      = scrna.markers_sig$pct.2,
                                  "FlyBase_ID" = paste0("<a  href='https://flybase.org/reports/",scrna.markers_sig$geneID,"'>", scrna.markers_sig$geneID,"</a>"))

dtable <- DT::datatable(final_df,class = 'cell-border stripe',rownames=F,filter='top',
              editable = TRUE, extensions = 'Buttons', options = list(
                dom = 'Bfrtip',
                buttons = c('copy','csv','excel','pdf','print')
              ),escape = FALSE)
}
saveWidget(dtable,"FindMarkers.html")
            
if (species == "Homo Sapiens"){
# Create Interactive Table
tibble::as_tibble(unique(scrna.markers_sig)) %>% dplyr::arrange(p_val_adj) -> scrna.markers_sig
scrna.markers_sig = scrna.markers_sig %>% dplyr::select(cluster, geneSymbol, p_val, avg_log2FC,
                                             p_val_adj,pct.1,pct.2,geneID)
final_df <- data.frame(           "cluster"    = scrna.markers_sig$cluster,
                                  "gene"       = scrna.markers_sig$geneSymbol,
                                  "p_val"      = scrna.markers_sig$p_val,
                                  "avg_log2FC" = scrna.markers_sig$avg_log2FC,
                                  "p_val_adj"  = scrna.markers_sig$p_val_adj,
                                  "pct.1"      = scrna.markers_sig$pct.1,
                                  "pct.2"      = scrna.markers_sig$pct.2,
                                  "GeneCards" = paste0("<a  href='https://www.genecards.org/Search/Keyword?queryString=",scrna.markers_sig$geneID,"'>", scrna.markers_sig$geneID,"</a>"))

dtable <- DT::datatable(final_df,class = 'cell-border stripe',rownames=F,filter='top',
              editable = TRUE, extensions = 'Buttons', options = list(
                dom = 'Bfrtip',
                buttons = c('copy','csv','excel','pdf','print')
              ),escape = FALSE)
}
saveWidget(dtable,"FindMarkers.html")

# Feature plot
pdf(paste0(out_dir, "combined.top25markers.pdf"))
ggp = list()
for (marker in topN){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker)
    print(ggp[[marker]])
}
dev.off()


#DoHeatmap(scrna.combined, features = topN$geneID, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend() + scale_y_discrete(breaks=topN$geneID,
#        labels=geneTable$geneSymbol[match(topN$geneID, geneTable$geneID)])
#ggsave(paste0(out_dir, "top25markers.heatmap.integrated.geneSymbol.pdf"), width = 8, height = 12)
#  

remove_markers <- setdiff(markers,row.names(scrna.combined))
markers <- markers[!markers%in%remove_markers]
markers <- unique(markers)
pdf(paste0(out_dir, "combined.markers.geneSymbol.pdf"))
ggp = list()
for (marker in markers){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker) + ggtitle(marker)
    print(ggp[[marker]])
}
dev.off()
  
pdf(paste0(out_dir, "combined.dotplot.geneSymbol.pdf"), width = 30, height = 10)
p1 <- DotPlot_scCustom(scrna.combined, features = markers, x_lab_rotate = TRUE, colors_use = "blue") + scale_x_discrete(breaks= markers)
print(p1)
dev.off()
      
# developer, clustered
pdf(paste0(out_dir, "combined.dotplot.clustered.geneSymbol.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = markers, plot_km_elbow = F, new_row_labels = markers)
print(p1)
dev.off()
