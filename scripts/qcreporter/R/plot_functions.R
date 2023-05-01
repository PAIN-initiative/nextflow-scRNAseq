#' QC Stacked Barplot with Faceting
#'
#' The metric used for category_x will generate bars as columns on the x-axis. The bar will be split vertically based on category_y.
#' Facet wrapping will be performed on supplied variables based on the facet_formula.
#'
#' @param meta A data.frame containing metadata
#' @param category_x A character object specifying the metadata to use for splitting in the x-direction
#' @param name_x A character object specifying a name to display on the x-axis
#' @param category_y  A character object specifying the metadata to use for color groups
#' @param category_name A character object specifying a name to display for the colors
#' @param colorset_y A colorset to use as fills for category_y. Currently supported: "rainbow" or "varibow". Default is "varibow"
#' @param name_y 	A character object specifying a name for the y-axis.
#' @param as_fraction A logical object specifying whether or not to display the stacked bars as fractions of the total count for each category_x. Default is FALSE.
#' @param facet_formula A formula object for faceting based on variables in the meta data frame. For example, \code{formula("~pool_id")} will facet wrap by a variable called pool_id in meta.
#' @param ... Additional arguments passed to \code{ggplot2::facet_wrap()}
#' @return A ggplot2 plot object
#' @import data.table
#' @import ggplot2
#' @export
qc_stacked_barplot_facet <- function (meta,
                                      category_x = "batch_id",
                                      name_x = "Batch ID",
                                      category_y = "well_id",
                                      category_name = "Well ID",
                                      colorset_y = "varibow",
                                      name_y = "N Cells",
                                      as_fraction = FALSE,
                                      facet_formula = NULL, ...) {
  assertthat::assert_that(sum(class(meta) %in% c("data.frame",
                                                 "data.table")) > 0)
  assertthat::assert_that(class(category_x) == "character")
  assertthat::assert_that(length(category_x) == 1)
  assertthat::assert_that(category_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(category_y) == "character")
  assertthat::assert_that(length(category_y) == 1)
  assertthat::assert_that(category_y %in% names(meta))
  assertthat::assert_that(class(category_name) == "character")
  assertthat::assert_that(length(category_name) == 1)
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(colorset_y) == "character")
  assertthat::assert_that(length(colorset_y) == 1)
  assertthat::assert_that(colorset_y %in% c("rainbow", "varibow"))
  assertthat::assert_that(class(as_fraction) == "logical")
  assertthat::assert_that(is.null(facet_formula) || rlang::is_formula(facet_formula))
  meta <- data.table::as.data.table(meta)

  if(!is.null(facet_formula)){
    formula_cols <- as.character(as.list(facet_formula))
    f_cols <- setdiff(formula_cols, "`~`")
    count_table <- meta[, .(n_cells = nrow(.SD)), by = mget(c(category_x,
                                                              category_y, f_cols))]
  } else {
    count_table <- meta[, .(n_cells = nrow(.SD)), by = mget(c(category_x,
                                                              category_y))]
  }
  plot_xpos <- data.frame(unique(count_table[[category_x]]))
  names(plot_xpos) <- category_x
  plot_xpos <- plot_xpos[order(plot_xpos[[category_x]]), ,
                         drop = FALSE]
  plot_xpos$xpos <- 1:nrow(plot_xpos)
  count_table <- count_table[plot_xpos, on = category_x]
  plot_fills <- data.frame(unique(count_table[[category_y]]))
  names(plot_fills) <- category_y
  if (colorset_y == "rainbow") {
    set.seed(3030)
    plot_fills$fill <- sample(grDevices::rainbow(nrow(plot_fills)),
                              nrow(plot_fills))
  }
  else if (colorset_y == "varibow") {
    set.seed(3030)
    plot_fills$fill <- sample(H5MANIPULATOR::varibow(nrow(plot_fills)),
                              nrow(plot_fills))
  }
  plot_fills <- plot_fills[order(plot_fills[[category_y]]),
  ]
  count_table <- count_table[plot_fills, on = category_y]
  count_table <- count_table[order(get(category_y), decreasing = TRUE)]
  if (as_fraction) {
    count_table <- count_table[, `:=`(ymax, cumsum(n_cells)/sum(n_cells)),
                               by = list(get(category_x))]
    count_table <- count_table[, `:=`(ymin, shift(ymax,
                                                  fill = 0, type = "lag")), by = list(get(category_x))]
  }
  p <- ggplot() + geom_rect(data = count_table,
                                              aes(xmin = xpos - 0.4, xmax = xpos + 0.4, ymin = ymin,
                                                           ymax = ymax, fill = fill)) + scale_fill_identity(category_name,
                                                                                                                     breaks = plot_fills$fill, labels = plot_fills[[category_y]],
                                                                                                                     guide = "legend") + scale_x_continuous(name_x,
                                                                                                                                                                     breaks = plot_xpos$xpos, labels = plot_xpos[[category_x]]) +
    scale_y_continuous(name_y) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                                       hjust = 1, vjust = 0.3))
  if(!is.null(facet_formula)){
    p <- p +
      facet_wrap(facet_formula, ...)
  }
  p
}

#' Generate a baseline-aligned barplot for two categorical metrics with faceting
#'
#' The metric used for category_x will generate bars as columns on the x-axis.
#' The bar will be split vertically based on category_y. Each group in category_y
#' will be aligned to make them easier to compare. Facet wrapping will be
#' performed on supplied variables based on the facet_formula. By specifying "stat"
#' the function can accept either metadata that has been pre-aggregated or raw values
#' to be counted.
#'
#' @param meta A data.frame containing metadata
#' @param category_x A character object specifying the metadata to use for grouping on the x-axis
#' @param name_x A character object specifying a name to display on the x-axis
#' @param category_y A character object specifying the metadata to use for splitting in the y-direction
#' @param category_name A character object specifying a name to display for the colors
#' @param colorset_y A colorset to use as fills for category_y. Currently supported: "rainbow" or "varibow". Default is "varibow"
#' @param name_y A character object specifying a name for the y-axis.
#' @param padding A numeric object specifying the fraction of the total vertical space to use for separating category_y groups. Default is 0.2.
#' @param stat A character object specifying the stat method for the barplot. Default is "count", also supports "identity".
#' @param variable_y_identity A character object specifying a the values to plot on y axis when stat="identity"
#' @param facet_formula A formula object for faceting based on variables in the meta data frame. For example, \code{formula("~pool_id")} will facet wrap by a variable called pool_id in meta.
#' @param ... Additional arguments passed to \code{ggplot2::facet_wrap()}
#'
#' @return a ggplot2 plot object
#' @export
#' @examples
#' set.seed(3)
#' test_data_1 <- data.frame(batch_id=sample(c("B001","B002"),1000, replace = TRUE),
#'                         well_id = sample(c("W1","W2","W3"),1000, replace = TRUE),
#'                         sample_names= 1:1000)
#' qc_aligned_barplot_facet(test_data_1,
#'                          category_x = "batch_id",
#'                          name_x = "Batch ID",
#'                          category_y = "well_id",
#'                          category_name = "Well ID",
#'                          colorset_y = "varibow",
#'                          name_y = "N Cells",
#'                          padding = 0.2)
#'
#' test_data <- rbind(data.frame(row=1:1000,
#'                               pool_id="P1",
#'                               well_id = sample(c("W1","W2","W3"),1000, replace = TRUE),
#'                               sample_id= sample(1:3, 1000, replace = TRUE)),
#'                    data.frame(row=1001:2000,
#'                               pool_id="P2",
#'                               well_id = sample(c("W4","W5","W6"),1000, replace = TRUE),
#'                               sample_id= sample(4:6, 1000, replace = TRUE)))
#' qc_aligned_barplot_facet(test_data,
#'                          category_x = "well_id",
#'                          name_x = "Well ID",
#'                          category_y = "sample_id",
#'                          category_name = "Sample ID",
#'                          colorset_y = "varibow",
#'                          name_y = "N Cells",
#'                          padding = 0.2,
#'                          facet_formula = formula("~pool_id"),
#'                          scales = "free")
#'
#' \dontrun{  # uses tidyr
#' test_data_3 <- rbind(tidyr::expand_grid(pool_id = "P1",
#'                                sample_id = c("S1","S2","S3"),
#'                                well_id = c("W1","W2","W3")),
#'                    tidyr::expand_grid(pool_id = "P2",
#'                                sample_id = c("S4","S5","S6"),
#'                                well_id = c("W4","W5","W6"))
#'                  )
#' test_data_3$n_cells <- round(rnorm(nrow(test_data_2), 1000, 100))
#' qc_aligned_barplot_facet(test_data_3,
#'     category_x = "well_id",
#'     name_x = "Well ID",
#'     category_y = "sample_id",
#'     category_name = "Sample ID",
#'     colorset_y = "varibow",
#'     name_y = "N Cells",
#'     padding = 0.2,
#'     stat= "identity",
#'     variable_y_identity = "n_cells",
#'     facet_formula = formula("~pool_id"),
#'     scales = "free")
#' }
qc_aligned_barplot_facet <- function(meta,
                                     category_x = "batch_id",
                                     name_x = "Batch ID",
                                     category_y = "well_id",
                                     category_name = "Well ID",
                                     colorset_y = "varibow",
                                     name_y = "N Cells",
                                     padding = 0.2,
                                     stat = "count",
                                     variable_y_identity = NULL,
                                     facet_formula = NULL,
                                     ...) {

  assertthat::assert_that(sum(class(meta) %in% c("data.frame","data.table")) > 0)
  assertthat::assert_that(class(category_x) == "character")
  assertthat::assert_that(length(category_x) == 1)
  assertthat::assert_that(category_x %in% names(meta))
  assertthat::assert_that(class(name_x) == "character")
  assertthat::assert_that(length(name_x) == 1)
  assertthat::assert_that(class(category_y) == "character")
  assertthat::assert_that(length(category_y) == 1)
  assertthat::assert_that(category_y %in% names(meta))
  assertthat::assert_that(class(category_name) == "character")
  assertthat::assert_that(length(category_name) == 1)
  assertthat::assert_that(class(name_y) == "character")
  assertthat::assert_that(length(name_y) == 1)
  assertthat::assert_that(class(colorset_y) == "character")
  assertthat::assert_that(length(colorset_y) == 1)
  assertthat::assert_that(colorset_y %in% c("rainbow","varibow"))
  assertthat::assert_that(class(padding) == "numeric")
  assertthat::assert_that(length(padding) == 1)
  assertthat::assert_that(padding < 1)
  assertthat::assert_that(length(stat) == 1)
  assertthat::assert_that(class(stat) == "character")
  assertthat::assert_that(stat %in% c("count","identity"), msg = "parameter stat must be 'count' or 'identity'")
  assertthat::assert_that(ifelse(stat == "identity", length(variable_y_identity) == 1,TRUE),
                          msg = "If stat is 'identity', variable_y_idnetity must be supplied")
  assertthat::assert_that(ifelse(stat == "identity", class(variable_y_identity) == "character",TRUE))
  assertthat::assert_that(ifelse(stat == "identity", variable_y_identity %in% names(meta),TRUE))
  assertthat::assert_that(is.null(facet_formula) || class(facet_formula) == "formula")

  meta <- as.data.table(meta)
  if(!is.null(facet_formula)){
    f_list <- as.character(as.list(facet_formula))
    f_cols <- setdiff(f_list, c("`~`","[.]","+"))
    f_cols <- trimws(unlist(strsplit(f_cols, split = "\\+")))
    assertthat::assert_that(all(f_cols %in% names(meta)), msg = "Columns in facet formula must be present in input meta object")

    if(stat == "count"){
      count_table <- meta[, .(counts = nrow(.SD)), by = mget(c(category_x, category_y, f_cols))]
    } else if (stat == "identity"){
      count_table <- meta[, mget(c(category_x,category_y, f_cols, variable_y_identity))]
      count_table[, counts := get(variable_y_identity)]
    }
  } else {
    if(stat == "count"){
      count_table <- meta[, .(counts = nrow(.SD)), by = mget(c(category_x, category_y))]
    } else if (stat == "identity"){
      count_table <- meta[, mget(c(category_x, category_y, variable_y_identity))]
      count_table[, counts := get(variable_y_identity)]
    }
  }

  plot_xpos <- data.frame(unique(count_table[[category_x]]))
  names(plot_xpos) <- category_x
  plot_xpos <- plot_xpos[order(plot_xpos[[category_x]]),,drop = FALSE]
  plot_xpos$xpos <- 1:nrow(plot_xpos)

  count_table <- count_table[plot_xpos, on = category_x]

  plot_fills <- data.frame(unique(count_table[[category_y]]))
  names(plot_fills) <- category_y
  if(colorset_y == "rainbow") {
    set.seed(3030)
    plot_fills$fill <- sample(grDevices::rainbow(nrow(plot_fills)), nrow(plot_fills))
  } else if(colorset_y == "varibow") {
    set.seed(3030)
    plot_fills$fill <- sample(H5MANIPULATOR::varibow(nrow(plot_fills)), nrow(plot_fills))
  }
  plot_fills <- plot_fills[order(plot_fills[[category_y]]),]
  count_table <- count_table[plot_fills, on = category_y]

  group_maxes <- count_table[, .(group_max = max(counts)), by = list(get(category_y))]
  names(group_maxes)[1] <- category_y
  group_maxes <- group_maxes[order(get(category_y), decreasing = TRUE)]
  group_maxes <- group_maxes[, cum_max := cumsum(group_max)]
  group_maxes <- group_maxes[, group_center := cum_max - group_max / 2]
  group_maxes <- group_maxes[, padded_center := group_center + (max(cum_max) * (padding/nrow(group_maxes))) * (1:nrow(group_maxes) - 1)]
  group_maxes <- group_maxes[, padded_base := padded_center - group_max/2]
  group_maxes <- group_maxes[, padded_top := padded_center + group_max/2]

  count_table <- count_table[group_maxes, on = category_y]

  count_table <- count_table[order(get(category_y), decreasing = TRUE)]
  count_table <- count_table[, ymax := cumsum(counts), by = list(get(category_x))]
  count_table <- count_table[, ymin := shift(ymax, fill = 0, type = "lag"), by = list(get(category_x))]

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = count_table,
                       ggplot2::aes(xmin = xpos - 0.4,
                                    xmax = xpos + 0.4,
                                    ymin = padded_base,
                                    ymax = padded_base + counts,
                                    fill = fill)) +
    ggplot2::geom_hline(data = count_table,
                        ggplot2::aes(yintercept = padded_base)) +
    ggplot2::geom_hline(data = count_table,
                        ggplot2::aes(yintercept = padded_top),
                        linetype = "dashed") +
    ggplot2::scale_fill_identity(category_name,
                                 breaks = plot_fills$fill,
                                 labels = plot_fills[[category_y]],
                                 guide = "legend") +
    ggplot2::scale_x_continuous(name_x,
                                breaks = plot_xpos$xpos,
                                labels = plot_xpos[[category_x]]) +
    ggplot2::scale_y_continuous(name_y,
                                breaks = c(group_maxes$padded_base,
                                           group_maxes$padded_top),
                                labels = c(rep("", nrow(group_maxes)),
                                           group_maxes$group_max),
                                expand = ggplot2::expand_scale(c(0, 0.02))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.3))

  if(!is.null(facet_formula)){
    p <- p +
      facet_wrap(facet_formula, ...)
  }

  p
}


#' Add Axis Spacing for Faceted Plotly
#'
#' Note: This does not work as expected in Rmarkdown
#' Add lines of whitespace to x-axis or y-axis title to shift title down or left.
#'
#'
#' @param plotly_obj An object generated by ggplotly
#' @param axis String value of "x" or "y". Which axis title to add spacing to.
#' If need to add spacing to both, make separate calls to `add_axis_title_spacing_plotly()`
#' @param n_lines Numeric value. Number of lines worth of spacing to add.
#' @return A the original plotly object with updated axis name
#' @examples
#' \dontrun{  # THIS FUNCTION NOT CURRENTLY EXPORTED
#' library(ggplot2)
#' library(plotly)
#' library(dplyr)
#' set.seed(1)
#' x <- sapply(1:10,function(x){paste(LETTERS[sample(1:26,10,replace = TRUE)],collapse="")})
#' y <- sapply(1:20,function(x){paste(LETTERS[sample(1:26,10,replace = TRUE)],collapse="")})
#' z <- rep(c("A","B"), times = 10)
#' df <- data.frame(x, y, z)
#' g <- ggplot(df, aes(x, y))+
#'   geom_point()+
#'   facet_wrap(~z, nrow =1) +
#'   xlab("X TITLE") +
#'   ylab("Y TITLE") +
#'   theme(axis.text.x = element_text(angle = 90))
#' g
#' gp <- ggplotly(g)
#' gp
#' gp %>%
#'   add_axis_title_spacing_plotly("x", 3) %>%
#'   add_axis_title_spacing_plotly("y", 3)
#' }
#'
add_axis_title_spacing_plotly <- function(plotly_obj, axis, n_lines){
  assertthat::assert_that(is.character(axis))
  axis <- tolower(axis)
  assertthat::assert_that(axis %in% c("x","y"))
  assertthat::assert_that(is.numeric(n_lines))

  axis_string <- paste0(axis,"axis")

  annotation_index <- (axis=="y") + 1
  orig_title <- plotly_obj$x$layout$annotations[[annotation_index]]$text
  if(is.null(orig_title)){
    stop(sprintf("No title detected in plotly object for axis %s", axis))
  }

  spacing_string <- paste("",rep("\n",n_lines), collapse = "")
  if(axis == "x"){
    newtitle <- paste0(spacing_string, orig_title)
  } else if (axis == "y"){
    newtitle <- paste0(orig_title, spacing_string)
  }

  plotly_obj$x$layout$annotations[[annotation_index]]$text <- newtitle

  return(plotly_obj)
}

#' Add Axis Spacing for Faceted Plotly
#'
#' Adjusts the axis titles in plotly objects by modifying the object itself,
#' since plotly's auto spacing does not do this appropriately for faceted plots
#'
#' @param plotly_obj An object generated by ggplotly
#' @param axis String value of "x" or "y". Which axis title to add spacing to.
#' If need to add spacing to both, make separate calls to `adjust_axis_title_spacing_plotly()`
#' @param adjustment Numeric value. Vertical (for x-axis) or horizontal (for y-axis) adjustment for axis title
#' @return A the original plotly object with updated axis name
#' @export
#' @examples
#' library(ggplot2)
#' library(plotly)
#' library(dplyr)
#' set.seed(1)
#' x <- sapply(1:10,function(x){paste(LETTERS[sample(1:26,10,replace = TRUE)],collapse="")})
#' y <- sapply(1:20,function(x){paste(LETTERS[sample(1:26,10,replace = TRUE)],collapse="")})
#' z <- rep(c("A","B"), times = 10)
#' df <- data.frame(x, y, z)
#' g <- ggplot(df, aes(x, y))+
#'   geom_point()+
#'   facet_wrap(~z, nrow =1) +
#'   xlab("X TITLE") +
#'   ylab("Y TITLE") +
#'   theme(axis.text.x = element_text(angle = 90))
#' g
#' gp <- plotly::ggplotly(g)
#' gp
#' gp %>%
#'   adjust_axis_title_spacing_plotly("x", 0.05) %>%
#'   adjust_axis_title_spacing_plotly("y", 0.05)
#'
adjust_axis_title_spacing_plotly <- function(plotly_obj, axis, adjustment){
  assertthat::assert_that(is.character(axis))
  axis <- tolower(axis)
  assertthat::assert_that(axis %in% c("x","y"))
  assertthat::assert_that(is.numeric(adjustment))

  annotation_index <- (axis=="y") + 1

  if(axis == "x"){
    plotly_obj$x$layout$annotations[[annotation_index]]$y <- adjustment
  } else if (axis == "y"){
    plotly_obj$x$layout$annotations[[annotation_index]]$x <- adjustment
  }

  return(plotly_obj)
}

#' Plot UMAP
#'
#' Standardized UMAP plotting with color customization for grouping variables
#'
#' @param df A data frame containing umap data and metadata for plotting
#' @param x_col Character value. Name of data column to plot on the x-axis
#' @param x_lab Character value. Name to display on the x-axis
#' @param y_col CHaracter value. Name of data column to plot on the y-axis
#' @param y_lab Character value. Name to display on the y-axis
#' @param title Character value. Title of plot.
#' @param point_size Numeric value. Point size to use for scatter plot
#' @param color_col Character value, Name of data column used to color point by groups or value
#' @param scale_color_fun A function that returns color scaling for ggplot2
#' @param ... Additional arguments passed to \code{scale_color_fun}
#' @return A ggplot2 plot object
#' @import ggplot2
#' @export
plot_umap_report <- function(df, x_col, x_lab, y_col, y_lab, title, point_size, color_col, scale_color_fun,...){
  g <- ggplot(df, aes_string(x_col, y_col)) +
    geom_point(alpha = 1, size = point_size, aes_string(color = color_col)) +
    ggtitle(title)+
    xlab(x_lab) +
    ylab(y_lab) +
    scale_color_fun(...) +
    theme_bw() +
    theme(aspect.ratio = 1/1,
          text = element_text(size = 20))
  return(g)

}

#' Pretty Color-Gradient-N Scaling across Data
#'
#' Creates a standard rainbow color scaling for set of values. Ensures default range
#' of interest will be colored consistently, with any observed values higher than
#' the default range colored a separate color. By default, expected range is blue to
#' red with high values as dark red. Default ranges were roughly based on genes per cell
#'
#' @param max_value Character value. Highest data value to extend color scale to.
#' @param colors Colors to use in  \code{scale_color_gradientn()}
#' @param value_breaks Values associated with each color, in order
#' @param highest_col Color to scale out to for observed values higher than
#' largest \code{value_breaks} value
#' @return A custom function of \code{scale_color_gradientn()} with colors, values, and breaks defined.
#' The function will accept additional arguments to \code{scale_color_gradientn()}
#' @import ggplot2
#' @export
#' @examples
#' library(ggplot2)
#' set.seed(3)
#' high_genes <- data.frame(x = runif(100, 0, 10), y = 1:100, n_genes = rnorm(100, 3500, 1000))
#' my_color_fun_high <- scale_color_manual_rainbow(max(high_genes$n_genes))
#' ggplot(high_genes, aes(x, y, color = n_genes))+ geom_point() + my_color_fun_high()
#'
#' super_high_genes <- data.frame(x = runif(100, 0, 10), y = 1:100, n_genes = runif(100, 3000, 9000))
#' my_color_fun_super_high <- scale_color_genes(max(super_high_genes$n_genes))
#' ggplot(super_high_genes, aes(x, y, color = n_genes))+ geom_point() + my_color_fun_super_high()
#'
#' low_genes <- data.frame(x = runif(100, 0, 10), y = 1:100, n_genes = runif(100, 0, 2900))
#' my_color_fun_low <- scale_color_genes(max(low_genes$n_genes))
#' ggplot(low_genes, aes(x, y, color = n_genes))+ geom_point() + my_color_fun_low()
scale_color_manual_rainbow <- function(max_value,
                                       colors = c("blue","deepskyblue","green3", "yellow","orange","red"),
                                       value_breaks = c(0,500, 1000, 2000, 3000, 4000),
                                       highest_col = "darkred"){
  function(...){
    if(max_value > max(value_breaks)){
      scale_color_gradientn(limits = c(0, max_value),
                            colours = c(colors, highest_col),
                            values = scales::rescale(c(value_breaks, max_value),
                                                     from = c(0, max_value)),
                            breaks = pretty(c(0, max_value), n = 4, min.n = 3), ...)
    } else {
      i_max_break <- which(value_breaks > max_value)[1] # smallest break greater than data
      scale_color_gradientn(limits = c(0, value_breaks[i_max_break]),
                            colours = colors[1:i_max_break],
                            values = scales::rescale(value_breaks[1:i_max_break],
                                                     from = c(0, value_breaks[i_max_break])),
                            breaks = pretty(c(0, value_breaks[i_max_break]), n = 4, min.n = 3), ...)
    }
  }
}

#' Default Color Gradient for Genes per Cell
#'
#' Uses \code{} with default values to create a gradient for genes per cell, with
#' scale adjusting to observed values
#'
#' @param max_value Character value. Highest genes per cell value in the dataset
#' @return A custom function of \code{scale_color_gradientn()} with colors, values,
#' and breaks defined. The function will accept additional arguments to \code{scale_color_gradientn()}
#' @export
scale_color_genes <- function(max_value){
  scale_color_manual_rainbow(max_value,
                             colors = c("blue","deepskyblue","green3", "yellow","orange","red"),
                             value_breaks = c(0, 500, 1000, 2000, 3000, 4000),
                             highest_col = "darkred")
}

#' Default Color Gradient for UMI per Cell
#'
#' Uses \code{} with default values to create a gradient for UMI per cell, with
#' scale adjusting to observed values
#'
#' @param max_value Character value. Highest UMI per cell value in the dataset
#' @return A custom function of \code{scale_color_gradientn()} with colors, values,
#' and breaks defined. The function will accept additional arguments to \code{scale_color_gradientn()}
#' @export
scale_color_umis <- function(max_value){
  scale_color_manual_rainbow(max_value,
                             colors = c("blue","deepskyblue","green3", "yellow","orange","red"),
                             value_breaks = c(0,1000, 3000, 5000, 7500, 10000),
                             highest_col = "darkred")
}


#' Default Color Gradient for Fraction Values
#'
#' Creates a default color gradient for fractional values using \code{scale_color_gradientn()}
#'
#' @param colors Colors for \code{scale_color_gradientn()}.
#' @param breaks Breaks for \code{scale_color_gradientn()}.
#' @param ... Additional arguments passed to \code{scale_color_gradientn()}
#' @return A custom function of \code{scale_color_gradientn()} for fractional values with colors, values,
#' and breaks defined. The function will accept additional arguments to \code{scale_color_gradientn()}
#' @export
scale_color_fct_mito <- function(colors = c("blue", "green3","yellow","red"),
                                 breaks = c(0,0.25,0.5,0.75,1), ...){
  ggplot2::scale_color_gradientn(limits = c(0, 1),
                       colors = colors,
                       breaks = breaks,...)
}

#' Seurat 3 Cell Palette
#'
#' Color palette with color assigned to each cell type in Seurat v3 labeling
#' reference (variable_pbmc_10k_v3)
#'
#' @return A data frame with columns "seurat3_pbmc_type" and "cell_color"
#' @export
seurat_3_cell_palette <- function(){
  data.frame(cell_labels = c("B cell progenitor", "CD14+ Monocytes", "CD16+ Monocytes",
                                   "CD4 Memory", "CD4 Naive","CD8 effector", "CD8 Naive",
                                   "Dendritic cell", "Double negative T cell", "NK cell",
                                   "pDC", "Platelets", "pre-B cell"),
             cell_color = c("#FF0000","#FF8C00","#FFEE00","#44FF00", "#00E1FF","#0000FF",
                            "#E546FA","#F598E5","#008A12", "#803CCF", "#967729", "#B1C4F0",
                            "#DCF0B1"),
             stringsAsFactors = FALSE)
}

#' Load Cell Label Color Palette
#'
#' For recognized labeling methods loads a cell-type specific color palette,
#' otherwise genreates a varibow palette based on number of cell types in data
#' @param cell_label_method Character value. Cell labeling method. Currently only
#' defined method is Seurat3.
#' @param cell_labels Character vector. All cell type labels in the dataset.
#' @return Dataframe of cell types and corresponding colors.
#' @export
#' @examples
#' get_cell_label_palette(cell_label_method = "Seurat3")
#' \dontrun{get_cell_label_palette(cell_label_method = "SpecialLabelMethod") ## generates intentional error}
#' get_cell_label_palette(cell_label_method = "SpecialLabelMethod",
#'                        cell_labels = c("Label1","Label2","Label3"))
get_cell_label_palette <- function(cell_label_method = "Seurat3", cell_labels=NULL){
  if(cell_label_method == "Seurat3"){
    seurat_3_cell_palette()
  } else {
    assertthat::assert_that(!is.null(cell_labels),
                            msg = "Unrecognized cell label method. Must supply cell_labels to generate a custom palette")
    unique_labels = unique(sort(cell_labels,na.last = TRUE))
    n_labels <- length(unique_labels)
    pal_colors <- H5MANIPULATOR::varibow(n_labels)
    df_pal <-   data.frame(cell_labels = unique_labels,
                           cell_color = pal_colors,
                           stringsAsFactors = FALSE)
    df_pal
  }
}

#' QC Pie Chart with Optional Faceting
#'
#' Create a pie chart with wedges specified by a grouping variable. Faceting may be performed
#' if facetting is provided. 
#'
#' The grouping variable used for category_y will generate the wedges of the bar plot. The supplied
#' dataframe can either be one observation per row or one row per category_y. If the former, "stat" should be 
#' "count", in which case counts of each category will be calculated. If the latter, "stat" should be
#' "identity", and a column containing the category counts should be supplied ("variable_y_identity"). 
#' Facet wrapping will be performed on supplied variables based on the facet_formula. 
#'
#' @param meta A data.frame containing metadata
#' @param category_y  A character value specifying the metadata column to use for color groups
#' @param category_name A character value specifying a name to display for the category colors
#' @param stat A character value, either "count" or "identity", depending on whether input data
#' has already been aggregated by the category_y grouping variable. Default is "count".
#' @param palette_df_group An optional data frame containing two columns defining the 
#' categorical data color palette. The first column contains label values for 
#' all levels of category_y and the second column contains the color for each 
#' label. If NULL (default) random colors will be generated via the colorset_y
#' parameter.
#' @param colorset_y A colorset to use as fills for category_y. Currently 
#' supported: "rainbow" or "varibow". Default is "varibow". Only used when
#' palette_df_group is NULL.
#' @param variable_y_identity Value to use as labels for the pie wedge. Currently supported
#' values are "p
#' @param plot_font_size Numeric value. General size of text on plot in font size. Default is 14.
#' @param label_size Size of text to label each pie wedge. Default is 3
#' @param text_repel_x Numeric value. X coordinate for plotting pie wedge labels 
#' using \code{ggrepel::geom_text_repel()}. Default is 1.5
#' @param nudge_x_label Numeric value. Distance to nudge label in X coordinate 
#' for plotting pie wedge labels using \code{ggrepel::geom_text_repel()}. Higher
#' values will increase line length. Default is 0.5.
#' @param facet_formula A formula object for faceting based on variables in the meta data frame. 
#' For example, \code{formula("~pool_id")} will facet wrap by a variable called pool_id in meta. 
#' Currently only actual column names (rather than variables) supported for this argument.
#' @param ... Additional arguments passed to \code{ggplot2::facet_wrap()}
#' @return A ggplot2 plot object
#' @import data.table
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples 
#' \dontrun{
#' set.seed(2021)
#' my_palette <- batchreporter::get_cell_label_palette()
#' cell_probs <- c(0.025, 0.20, 0.03, 0.192, 0.18, 0.14, 0.05,
#'                 0.005, 0.02, 0.1 ,0.001,0.008,0.049)
#' test_data_1 <- data.frame(cell_type = sample(my_palette[,1], 500, 
#'                                              prob = cell_probs, replace = TRUE),
#'                         group = sample(c("MT UMI>0.1", "MT UMI<0.1"), 500, 
#'                                        prob = c(0.1, 0.9), replace = TRUE),
#'                         sample_names= 1:500)
#' # Use a predefined palette
#' qc_piechart(test_data_1,
#'             category_y = "cell_type", 
#'             category_name = "Cell Type",
#'             stat = "count",
#'             palette_df_group = my_palette,
#'             plot_font_size = 14,
#'             label_size = 3,
#'             text_repel_x = 1.5, 
#'             nudge_x_label = 0.2,
#'             facet_formula = as.formula("~group"))
#' # Use default random palette
#' qc_piechart(test_data_1,
#'             category_y = "cell_type",
#'             stat = "count",
#'             label_size = 3,
#'             text_repel_x = 1.5, 
#'             nudge_x_label = 0.2,
#'             facet_formula = as.formula("~group"))
#' 
#' # No Facetting
#' qc_piechart(test_data_1,
#'             category_y = "cell_type",
#'             stat = "count",
#'             label_size = 3,
#'             text_repel_x = 1.5, 
#'             nudge_x_label = 0.2)
#'
#' # Stat = "identity" version, precalculated counts
#' test_data_2 <- data.frame(cell_type = rep(my_palette[,1], times = 2),
#'                          group = rep(c("MT UMI>0.1", "MT UMI<0.1"), each = nrow(my_palette)),
#'                          sample_names= 1:(2*nrow(my_palette)),
#'                          group_count = c(round(1000*cell_probs),round(4000*sample(cell_probs))))
#' test_data_2
#' qc_piechart(test_data_2,
#'             category_y = "cell_type",
#'             stat = "identity",
#'             variable_y_identity = "group_count",
#'             label_size = 3,
#'             text_repel_x = 1.5, 
#'             nudge_x_label = 0.2,
#'             facet_formula = as.formula("~group"))
#' } # end don't run


qc_piechart <- function(meta, 
                         category_y,
                         category_name = "Cell Type", 
                         stat = "count",
                         variable_y_identity = NULL,
                         palette_df_group = NULL, 
                         colorset_y = "varibow",
                         plot_font_size = 14, 
                         label_size = 3,
                         text_repel_x = 1.5,
                         nudge_x_label = 0.5,
                         facet_formula = NULL, ...) {
  
  assertthat::assert_that(sum(class(meta) %in% c("data.frame",
                                                 "data.table")) > 0)
  assertthat::assert_that(class(category_y) == "character")
  assertthat::assert_that(length(category_y) == 1)
  assertthat::assert_that(category_y %in% names(meta))
  assertthat::assert_that(class(category_name) == "character")
  assertthat::assert_that(length(category_name) == 1)
  assertthat::assert_that(ifelse(!is.null(palette_df_group),
                                 all(unlist(meta[,category_y]) %in% palette_df_group[,1]),
                                 TRUE), 
                          msg = "Supplied color palette must contain all levels of category variable.")
  assertthat::assert_that(class(colorset_y) == "character")
  assertthat::assert_that(length(colorset_y) == 1)
  assertthat::assert_that(colorset_y %in% c("rainbow", "varibow"))
  assertthat::assert_that(length(stat) == 1)
  assertthat::assert_that(class(stat) == "character")
  assertthat::assert_that(stat %in% c("count", "identity"), msg = "parameter stat must be 'count' or 'identity'")
  assertthat::assert_that(ifelse(stat == "identity", length(variable_y_identity) ==
                                   1, TRUE), msg = "If stat is 'identity', variable_y_identity must be supplied")
  assertthat::assert_that(ifelse(stat == "identity", class(variable_y_identity) ==
                                   "character", TRUE))
  assertthat::assert_that(ifelse(stat == "identity", variable_y_identity %in% names(meta), TRUE))
  assertthat::assert_that(is.null(facet_formula) || 
                            class(facet_formula) == "formula")
  meta <- as.data.table(meta)
  
  # Format and calculate count table using relevant grouping variables and values
  if (!is.null(facet_formula)) {
    f_list <- as.character(as.list(facet_formula))
    f_cols <- setdiff(f_list, c("`~`", "[.]", "+"))
    f_cols <- trimws(unlist(strsplit(f_cols, split = "\\+")))
    assertthat::assert_that(all(f_cols %in% names(meta)), 
                            msg = "Columns in facet formula must be present in input meta object")
    if (stat == "count") {
      count_table <- meta[, .(counts = nrow(.SD)), by = mget(c(category_y, f_cols))]
    }
    else if (stat == "identity") {
      count_table <- meta[, mget(c(category_y, f_cols, variable_y_identity))]
      count_table[, `:=`(counts, get(variable_y_identity))]
    }
    count_table[, `:=`(pct = counts/sum(counts)) , by = mget(c(f_cols))]
    count_table <- count_table[base::order(count_table[, ..category_y], decreasing = TRUE),]
    count_table[, `:=`(ypos =  cumsum(pct)- 0.5*pct) , by = mget(c(f_cols))]
  }
  else {
    if (stat == "count") {
      count_table <- meta[, .(counts = nrow(.SD)), by = mget(category_y)]
    }
    else if (stat == "identity") {
      count_table <- meta[, mget(c(category_y, variable_y_identity))]
      count_table[, `:=`(counts, get(variable_y_identity))]
    }
    count_table[, `:=`(pct = counts/sum(counts))]
    count_table <- count_table[base::order(count_table[, ..category_y], decreasing = TRUE),]
    count_table[, `:=`(ypos =  cumsum(pct)- 0.5*pct) ]
    
  }
  
  # Color palette
  if(!is.null(palette_df_group)){
    plot_fills <- palette_df_group
    names(plot_fills) <- c(category_y, "fill")
  } else {
    plot_fills <- data.frame(unique(count_table[[category_y]]))
    names(plot_fills) <- category_y
    if (colorset_y == "rainbow") {
      set.seed(3030)
      plot_fills$fill <- sample(grDevices::rainbow(nrow(plot_fills)), 
                                nrow(plot_fills))
    }
    else if (colorset_y == "varibow") {
      set.seed(3030)
      plot_fills$fill <- sample(H5MANIPULATOR::varibow(nrow(plot_fills)), 
                                nrow(plot_fills))
    }
    plot_fills <- plot_fills[order(plot_fills[[category_y]]), ]
  }
  
  # Plot
  p <- ggplot(count_table, aes(x = 1, y=pct, fill=!!rlang::parse_expr(category_y))) +
    geom_bar(stat="identity", width = 1, color="white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(name = category_name, values = plot_fills$fill, breaks = plot_fills[,category_y]) +
    ggrepel::geom_text_repel(aes(label = !!rlang::parse_expr(category_y), x=text_repel_x, y = ypos), 
                    color = "black", size = label_size, segment.size = 0.4, nudge_x = nudge_x_label,
                    show.legend = TRUE) +
    theme_void() + 
    theme(text = element_text(size = plot_font_size), 
          strip.text = element_text(face = "bold"))
  
  if (!is.null(facet_formula)) {
    p <- p + facet_wrap(facet_formula, ...)
  }
  
  return(p)
}
