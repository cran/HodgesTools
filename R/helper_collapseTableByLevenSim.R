#' helper_collapseTableByLevenSim
#'
#' Reads in a table and value for Levenshtein threshold and returns a table collapsed by threshold (highest p-value for each group)
#'
#' @param inputTable dataframe. HOMER output table modified in the parent script--ready for filtering by Levenshtein similarity.
#' @param levenSimThresholdVal float. Value for thresholding TFs. For groups of TFs with similar consensus sequences, the TF with the lowest p-value by HOMER will be retained.
#'
#' @author Tim Scott
#'
#' @return tibble
#'


#### collapseTableByLevenSim
# Takes: A HOMER TF Known Results Table
# Returns: The same table, but collapsed by Levenshtein similarity
#
# Test variables:
# inputTable=inputListFiles_qValFilt$Liver
# levenSimThresholdVal=.9
# Assumes the table is a basic HOMER output table minus crap stock headers
helper_collapseTableByLevenSim <- function(inputTable,levenSimThresholdVal){
  # Get a 2-col DF: consensus and max Levenshtein Sim value of consensus row and all those above
  LevenSimCollapsedMotifDF <- helper_getMaxLevenSimCol(inputTable$Consensus)
  # Use second column to filter out consensuses above threshold (those which are represented by an earlier family consensus)
  LevelCollapsedMotifVector_underThreshold <- LevenSimCollapsedMotifDF[LevenSimCollapsedMotifDF$LevenSim<levenSimThresholdVal,1]
  # Return the input table filtered by the vectorlist of levenSim collapsed consensuses
  outTable <- dplyr::filter(inputTable, .data$Consensus %in% LevelCollapsedMotifVector_underThreshold)
  return(outTable)
}
