#' makeBigTableFromListofTables
#'
#' Reads in a list of tables and return list of tables with percent Fold Change (enrichment)
#'
#' @param inputListOfTables list of dataframe. List of HOMER TF knownResults Tables.
#'
#' @author Tim Scott
#'
#' @return list of tibble
#'

#### makeBigTableFromListofTables
# Takes: A basic HOMER TF known results table minus the header
# Returns: The same data frame but with cleaned percentage (char to num) and motif cols
helper_makeBigTableFromListOfStandardTables <- function(inputListOfTables){
  # Convert list to a single table for easier processing of all rows
  inputTables <- data.table::rbindlist(inputListOfTables, idcol="ID")
  # Change to factor to maintain order in plotting
  inputTables$ID <- factor(inputTables$ID, levels = unique(inputTables$ID))
  # Convert percentages and calculate fold change of %target / %background
  inputTables_withPercentCols <- inputTables %>%
    dplyr::mutate(percentTargetNum = (as.numeric(gsub("%", "", .data$percentTarget, fixed = TRUE))/100)) %>%
    dplyr::mutate(percentBackgroundNum = (as.numeric(gsub("%", "", .data$percentBackground, fixed = TRUE))/100)) %>%
    dplyr::mutate(percentFold = (.data$percentTargetNum/.data$percentBackgroundNum))
  # Convert the TF names to somthing more manageable
  inputTables_withPercentMotifNameCleaned <- inputTables_withPercentCols %>%
    dplyr::mutate(Motif = gsub("\\(.*", "", .data$MotifName))
  # Return
  return(inputTables_withPercentMotifNameCleaned)
}
