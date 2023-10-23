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
              default = NULL,
              help = "Input directory containing h5 and json files",
              metavar = "character"),
  make_option(opt_str = c("-z","--cellbender_dir"),
              type = "character",
              default = NULL,
              help = "cellbender_dir out directory",
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
              metavar = "character"),
  make_option(opt_str = c("-p","--alg"),
              type = "integer",
              default = 1,
              help = "Cluster Algorithm",
              metavar = "integer"))

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

if(is.null(args$experiment_id)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

if(!dir.exists(args$out_dir)) {
  dir.create(args$out_dir)
}

rmd_path <- file.path(args$out_dir,
                      paste0(args$experiment_id,
                             "_ngs_sample_qc_report.rmd"))

file.copy(system.file("rmarkdown/ngs_sample_qc_report.rmd", package = "qcreporter"),
          rmd_path,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_path,
  params = list(  experiment_id    = args$experiment_id,
                  cellbender_dir  = args$cellbender_dir,
                  outdir           = args$out_dir,
                  refdir           = args$refdir,
                  projectName      = args$experiment_id,
                  in_method_string = args$in_method,
                  in_dir           = args$in_dir,  
                  in_key           = args$in_key,  
                  out_dir          = args$out_dir,
                  species          = args$species,
                  resolution       = args$resolution,
                  percent_mito     = args$percent_mito,
                  percent_ribo     = args$percent_ribo,
                  filter_MALAT     = args$filter_MALAT,
                  filter_MITO      = args$filter_MITO,
                  filter_RIBO      = args$filter_RIBO,
                  alg              = args$alg),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_path)
