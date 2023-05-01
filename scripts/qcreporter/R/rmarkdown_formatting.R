# Rmarkdown formatting and QC utility functions

#' Knit a plot as a subchunk
#'
#' Knit a plot as a subchunk so its dimensions can be specified individually. Based
#' on code found here: http://michaeljw.com/blog/post/subchunkify/
#'
#' Allows individual plots within the same chunk to be knit as subchunks with
#' unique chunk options. Plots can be output in a loop with each plot using
#' different dimensions, ie dynamic dimensions based on number of x and/or y
#' category levels. Parent chunk should have chunk options 'results = "as-is"'
#' to ensure proper output. Note that this will create a "figures" directory in
#' the directory containing the Rmarkdown file containing the chunk plots. Ensure
#' that Rmarkdown  yaml has 'self_contained: true' in html document section (or equivalent)
#'
#' @param g The plot object
#' @param subchunk_name Character value. Unique name of Rmarkdown subchunk to be generated.
#' @param quiet_knit Logical value, default TRUE. Passed to \code{knitr::knit()}, should the subchunk
#' be knit "quietly" (no progress bar or messages)
#' @param chunk_opt_list Named list of chunk options for the subchunk. Can take any chunk
#' options available to a normal chunk.
#' @export
#' @examples
#' \dontrun{
#' # This will generate a file in 'figures' subdirectory of working directory
#' library(ggplot2)
#' g_example <- ggplot(data.frame(x=1:10, y = 1:10), aes(x, y)) + geom_point()
#' chunk_opt_l <- list(fig.height=10, fig.width=12, warning=TRUE)
#' make_subchunk(g_example, "test_chunk", chunk_opt_list = chunk_opt_l)
#' }
#'

make_subchunk <- function(g, subchunk_name, quiet_knit = TRUE, chunk_opt_list = list(fig.height=7, fig.width=5, warning = TRUE)) {
  if(is.null(subchunk_name)){
    subchunk_name <- paste0(as.numeric(Sys.time()), "_",)
  }

  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')

  # construct chunk option string
  if(!is.null(chunk_opt_list)){
    option_names <- names(chunk_opt_list)
    option_string_list <- sapply(1:length(chunk_opt_list), function(i){
      val <- chunk_opt_list[[i]]
      val_type <- class(val)
      quote_string <- ifelse(val_type=="character","'","")
      val_fmt <- paste0(quote_string, val, quote_string)
      paste(names(chunk_opt_list)[i], val_fmt, sep = "=")
    })
    option_string <- paste(c(", ",option_string_list), collapse = ", ")
  } else {
    option_string <- ""
  }

  # construct full chunk
  sub_chunk <- paste0("\n```{r ", subchunk_name, option_string, "}",
                      "\n(", g_deparsed, ")()",
                      "\n```")

  # knit chunk
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = quiet_knit))
}

#' Create a simple HTML table of key and value pairs
#'
#' Create a simple table of paired values (ie name and value).
#'
#' @param labels Vector of labels
#' @param values Vector of values, same length as labels
#' @param col_widths_px Numeric vector of length 2. Widths in pixels of each column,
#' @param fontsize Numeric value. Default 2.
#' @export
#' @return HTML of formatted table.
#' @examples
#' lab1 <- c("Batch", "N Samples", "N Pools")
#' val1 <- c("B000", 16, 2)
#' h_tab <- simple_html_table(lab1, val1)
#' htmltools::html_print(h_tab)
simple_html_table <- function(labels, values, col_widths_px = c(300, 800), fontsize = 2){
  assertthat::assert_that(length(col_widths_px) %in% c(1,2))
  assertthat::assert_that(length(labels) == length(values))
  assertthat::assert_that(all(col_widths_px > 100) & all(col_widths_px < 1024))
  # assertthat::assert_that(ifelse(is.null(col.names),TRUE, length(col.names) == 2))
  # assertthat::assert_that(is.null(col.names) | class(col.names) == "character")

  if(length(col_widths_px) == 1){
    col_widths_px <- c(col_widths_px, col_widths_px)
  }

  create_row <- function(value1, value2 , col_widths_px = c(50, 50)){
    sprintf('<tr>
                <td style="width:%spx">
                  %s
                </td>
                <td style="width:%spx">
                  %s
                </td>
            </tr>',
            col_widths_px[1], value1, col_widths_px[2], value2)
  }

  add_table_wrap <- function(all_row_string, fontsize = 2, col_widths_px= c(50,50)){
    sprintf('<font size="%s">
            <style>
              table {
                border: 1px solid black;
                border-collapse: collapse;
              }
              td {
                padding: 5px;
                border: 1px solid black;
                border-collapse: collapse;
              }
            </style>
            <table>
              <colgroup>
                <col span="1" style="width: %spx;">
                <col span="1" style="width: %spx;">
              </colgroup>
              <tbody>
                %s
              </tbody>
            </table>
            </font>', fontsize, col_widths_px[1], col_widths_px[2], all_row_string)
  }

  row_strings <- character()
  for (i in seq_along(labels)){
    row_strings[i] <- create_row(labels[i], values[i], col_widths_px = col_widths_px)
  }
  row_string_all <- paste(row_strings, collapse = "")

  out_table <- add_table_wrap(row_string_all,fontsize = fontsize, col_widths_px = col_widths_px )

  htmltools::HTML(out_table)
}


#' Format Flags in GT Table
#'
#' CURRENTLY NOT USED. Apply standard formatting to "Comments" column of a gt table
#'
#' For "Comments" column, will fill cells containing pattern "Warning" in light red
#' and will fill cells containing pattern "Fail" in dark red.
#'
#' @param x A gt table object. Must contain a column named "Comments" with flags "Warning" or "Fail".
#' @return A gt table object with formatted Comments cells.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(1)
#' vComments <- sample(c("Pass","Warning","Fail",NA), 10, replace = TRUE)
#' my_gt <- gt::gt(data.frame(X = 1:10, Comments = vComments))
#' my_gt
#' gt_fmt_comments(my_gt)
#' my_gt_nocommments <- gt::gt(data.frame(X = 1:10, CommentsX = vComments))
#' my_gt_nocommments
#' \dontrun{gt_fmt_comments(my_gt_nocommments)  # This intentionally generates an error}
gt_fmt_comments <- function(x) {
  assertthat::assert_that("Comments" %in% names(x[["_data"]]), msg = "Expect that gt contains flagging column named 'Comments'")

  x %>%
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "red" , alpha = 0.3)
      ),
      locations = gt::cells_body(
        columns = gt::vars(Comments),
        rows = grepl("Warning", .data$Comments))
    ) %>%
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "red" , alpha = 0.5)
      ),
      locations = gt::cells_body(
        columns = gt::vars(Comments),
        rows = grepl("Fail", .data$Comments))
    )
}


# Data summary functions

#' Get median and range
#'
#' Gets median and range of numeric vector, rounding to specified digits with option
#' to show the number of missing values. Median will always be calculated using
#' non-missing values.
#'
#' @param values Numeric vector
#' @param digits_round Number of digits to round to. Default 1.
#' @param comma_separate Logical, default TRUE. Should large values be displayed with comma separators
#' @param add_missing Logical, default TRUE. Should the number of missing values be displayed in brackets
#' @param verbose Logical, default TRUE. Should a message be printed upon successful calculation, including the number of missing values.
#' @return String value of formatted median and range. 'median (min-max) [missing]'
#' @import stats
#' @export
#' @examples
#' my_values <- c(1000:1010,NA)
#' get_median_range(my_values,
#'                  digits_round = 1,
#'                  comma_separate = TRUE,
#'                  add_missing = TRUE,
#'                  verbose = TRUE)
get_median_range <- function(values, digits_round = 1, comma_separate = TRUE, add_missing = TRUE, verbose = TRUE){
  assertthat::assert_that(mode(values) == "numeric" | all(is.na(values))) # can be any type of numeric, allow calculation if all missing
  assertthat::assert_that(length(values) > 0)

  fmt_num <- function(num){
    formatC(num, digits = digits_round, big.mark = ifelse(comma_separate,",",""), format = "f")
  }

  i_missing <- which(is.na(values))
  n_missing <- length(i_missing)
  has_missing <- n_missing > 0

  if(has_missing){
    values <- values[-i_missing]
  }

  med_string <- fmt_num(median(values))

  if(length(values) > 0){
    range_val <- fmt_num(range(values))
  } else {
    range_val <- c(NA,NA)
  }
  range_string <- paste(range_val, collapse = "-")
  med_range_string <- sprintf("%s (%s)", med_string, range_string)

  if(add_missing & has_missing){
    med_range_string <- sprintf("%s [%s]", med_range_string, n_missing)
  }

  if(verbose){
    cat(sprintf("Median and range calculated. Removed %s missing values", n_missing), sep = "\n")
  }

  return(med_range_string)
}


#' Get formatted range
#'
#' Gets range of numeric vector, rounding to specified digits with option
#' to show the number of missing values.
#'
#' @param values Numeric vector
#' @param digits_round Number of digits to round to. Default 1.
#' @param comma_separate Logical, default TRUE. Should large values be displayed with comma separators
#' @param add_missing Logical, default TRUE. Should the number of missing values be displayed in brackets
#' @param verbose Logical, default TRUE. Should a message be printed upon successful calculation, including the number of missing values.
#' @return String value of formatted median and range. 'median (min-max) [missing]'
#' @export
#' @examples
#' my_values <- c(1000:1010,NA)
#' get_range(my_values, digits_round = 1, comma_separate = TRUE, add_missing = TRUE, verbose = TRUE)
get_range <- function(values, digits_round = 1, comma_separate = TRUE, add_missing = TRUE, verbose = TRUE){
  assertthat::assert_that(mode(values) == "numeric" | all(is.na(values))) # can be any type of numeric, allow calculation if all missing
  assertthat::assert_that(length(values) > 0)

  fmt_num <- function(num){
    formatC(num, digits = digits_round, big.mark = ifelse(comma_separate,",",""), format = "f")
  }

  i_missing <- which(is.na(values))
  n_missing <- length(i_missing)
  has_missing <- n_missing > 0

  if(has_missing){
    values <- values[-i_missing]
  }

  if(length(values) > 0){
    range_val <- fmt_num(range(values))
  } else {
    range_val <- c(NA,NA)
  }
  range_string <- paste(range_val, collapse = "-")

  if(add_missing & has_missing){
    range_string <- sprintf("%s [%s]", range_string, n_missing)
  }

  if(verbose){
    cat(sprintf("Range calculated. Removed %s missing values", n_missing), sep = "\n")
  }

  return(range_string)
}

#' Determine passing specification
#'
#' This function is currently unused. Confirm functionality and modify if needed before use.
#' Flags numeric values above and/or below supplied upper and lower thresholds.
#'
#' @param values Numeric vector of values to evaluate for spec
#' @param lower_threshold Numeric value that represents the lower limit of acceptable values
#' @param upper_threshold Numeric value that represents the upper limit of acceptable values
#' @param digits_round The number of digits that values are to be rounded to before evaluation
#' @param pass_at_threshold Logical value. Whether or not values equal to threshold should pass, defaults to TRUE.
#' @param flag_values Vector of length 2 where the first value is the value to be returned for passing/non-flagged values
#' and the second value is the value to be returned for non-passing/flagged values.
#' @return Vector of length equal to input values. Output types defined by flag_values parameter
#' @export
#' @examples
#' myvalues <- 1:20
#' flags <- determine_passing_spec(myvalues,
#' lower_threshold = 5,
#' upper_threshold = 18,
#' pass_at_threshold = TRUE)
#' flags
#' flags2 <- determine_passing_spec(myvalues,
#' lower_threshold = 5,
#' upper_threshold = 18,
#' pass_at_threshold = FALSE)
#' flags2
determine_passing_spec <- function(values, lower_threshold = NULL, upper_threshold = NULL, digits_round = NULL, pass_at_threshold = TRUE, flag_values = c(FALSE,TRUE)){
  assertthat::assert_that(mode(values) == "numeric")
  assertthat::assert_that(length(values) >= 1)

  if (is.null(digits_round)){
    values_round <- values
  } else {
    values_round <- round(values, digits_round)
  }

  qc_flags <- rep(flag_values[1],length(values))

  # Flag values < or <= threshold
  if (!is.null(lower_threshold)){
    if (pass_at_threshold){
      i_low <- which(values_round < lower_threshold)
    } else {
      i_low <- which(values_round <= lower_threshold )
    }

    if (length(i_low) > 0){
      qc_flags[i_low] <- flag_values[2]
    }
  }

  # Flag values > or >= threshold
  if (!is.null(upper_threshold)){
    if (pass_at_threshold){
      i_high <- which(values_round > upper_threshold)
    } else {
      i_high <- which(values_round >= upper_threshold )
    }

    if (length(i_high) > 0){
      qc_flags[i_high] <- TRUE
    }
  }

  return(qc_flags)

}

markdown_colorize_html <- function(x, color) {
    sprintf("<span style='color: %s;'>%s</span>", color, x)
}

