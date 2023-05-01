#' Sample cells
#'
#' Samples cells on labels, retaining an equal number of each label type.
#'
#' @param cl Named character vector. Values are cell labels and names are cell barcodes.
#' @param sample.size Integer value. Number of cells of each cell type to sample
#' @param seed Integer value. Seed passed to \code{set.seed()} for reproducible sampling.
#' @return Named character vector. Values are the cell barcodes selected and names are the labels
#' @export
sample_cells <- function(cl,
                         sample.size,
                         seed = 3030) {
  cl.cells <- split(names(cl), cl)
  sampled.cells <- sapply(names(cl.cells), function(x) {
    cells <- cl.cells[[x]]
    if (sample.size >= length(cells)) {
      return(cells)
    }
    to.sample <- pmin(sample.size, length(cells))

    set.seed(seed)
    sample(cells, to.sample)
  }, simplify = FALSE)
  sampled.cells <- unlist(sampled.cells)
  return(sampled.cells)
}
