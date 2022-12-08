#' Read bed file
#'
#' Reads in a tab-delimited BED formatted file into R.
#'
#' First three columns of file must be the genomic coordinates of the regions (i.e. chr start end).
#'
#' read_bed will auto-detect BED3 and BED6 formats. It will also detect BED3+ and BED6+ formats assigning generic or user-defined col_names to the additional column(s).
#'
#'
#' @param file bed file
#' @param extra_col_names list of strings specifying extra column names
#' @param length boolean of whether to add length column
#' @param verbose boolean set to see function behavior
#' @author Tyler Hansen & Tim Scott
#' @examples
#' #load external data.
#' BED3 <- system.file(package = "HodgesTools", "extdata", "test_BED3.bed")
#' BED6 <- system.file(package = "HodgesTools", "extdata", "test_BED6.bed")
#' BED4 <- system.file(package = "HodgesTools", "extdata", "test_BED4.bed")
#' BED8 <- system.file(package = "HodgesTools", "extdata", "test_BED8.bed")
#'
#' # Read 3-column BED file.
#' read_bed(BED3)
#'
#' # Read 6-column BED file.
#' read_bed(BED6)
#'
#' # Read 3-column BED file and add length column.
#' read_bed(BED3, length = TRUE)
#'
#' # Read 3 column format BED file with additional fourth column. Add generic column names.
#' read_bed(BED4)
#'
#' # Read 3 column format BED file with additional fourth column. Specify additional column names.
#' read_bed(BED4, extra_col_names = c("fourthColumn"))
#'
#' # Read 6 column format BED file with additional columns. Specify additional column names.
#' read_bed(BED8, extra_col_names = c("seventhColumn", "eigthColumn"))
#'
#' @return tibble
#'
#' @export
read_bed <-
function(file, extra_col_names = c(), length = FALSE, verbose = TRUE) {
  ####


  options(readr.show_col_types = FALSE)

  #### Identify how many columns we are working with.
  firstLine <- suppressMessages(readr::read_tsv(file, n_max = 1, col_names = FALSE))
  nColFound <- length(firstLine)
  #### Stop if file is not BED format
  if (nColFound < 3) {
    stop("User-supplied file is fewer than three columns")
  }
  if (!(typeof(firstLine[[1]]) == "character" && (typeof(firstLine[[2]]) %in% c("double", "integer")) && (typeof(firstLine[[3]]) %in% c("double", "integer")))) {
    stop("User-supplied file is not in BED format.")
  }
  if (!grepl("^(c|C)hr([1][1-9]$|[0-9XYM]$|[2][1-2]$)", firstLine[[1]])) {
    stop("User-supplied file is not in BED format.")
  }

  ####
  #### Format column names for import
  # Set default column names
  baseColNames <- c("chr", "start", "end")
  isStandard6ColBed <- nColFound == 6 && typeof(firstLine[[5]]) %in% c("double", "integer") && firstLine[[6]] %in% c("+", "-", ".")

  # Perform checks to find input state
  if (nColFound == 3){ # (1) Check to see if this is a simple 3-col case
    msg <- "User-supplied file read as 3-column BED."
    cnamesToAssign <- baseColNames
  } else if (is.null(extra_col_names)) { # (2) Check to see if user supplied names
    # If no user input but >3-col and looks like standard 6-col BED, set BED6 names
    if (isStandard6ColBed) {
      msg <- "User-supplied file read as 6-column BED."
      cnamesToAssign <- c(baseColNames, c("name", "score", "strand"))
    } else { # Set additional default "col4","col5",[...] colnames
      msg <- "User-supplied file read as 3-colum BED with additional columns. Generic col_names added to 4th-Nth columns. User can specify col_names with extr_col_names argument."
      cnamesToAssign <- c(baseColNames, paste0("col", seq(1:nColFound)[4:nColFound]))
    }
  } else if (nColFound < 6) { #parse BED3+ vs BED6+
      if (length(extra_col_names) != (nColFound-3)) {
      # Check to make sure the amount of extra column names (extra_col_names) reflects the file.
      stop("User-supplied column names does not match column count.")
    } else { # Not BED3 and user input correct amount of names
      msg <- "User-supplied file read as 3-colum BED with additional columns. User-specified col_names added to 4th-Nth columns."
      cnamesToAssign <- c(baseColNames, extra_col_names)
    }
  } else if (nColFound > 6) {
    if (length(extra_col_names) != (nColFound-6)) {
      # Check to make sure the amount of extra column names (extra_col_names) reflects the file.
      stop("User-supplied column names does not match column count.")
    } else { # Not BED6 and user input correct amount of names
      msg <- "User-supplied file read as 6-colum BED with additional columns. User-specified col_names added to 7th-Nth columns."
      cnamesToAssign <- c(baseColNames, extra_col_names)
    }
  }

  #report message if verbose is TRUE.
  if (verbose) {
    message(msg)
  }

  ####
  #### Begin return
  # If user wants length, add length column. Else, return df.
  df <- suppressMessages(readr::read_tsv(file = file, col_names = cnamesToAssign))
  if (length) {
    df <- dplyr::mutate(df, length=(.data$end-.data$start))
  }
  # Return
  return(df)
}
