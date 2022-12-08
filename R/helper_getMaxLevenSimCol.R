#' helper_getMaxLevenSimCol
#'
#' Reads in a vector of motifs and returns a
#'
#' @param vectorOfMotifs vector of char. Vector of motifs to filter through.
#'
#' @author Tim Scott
#'
#' @return data.frame
#'

# This function is the meat for finding the similarity of one consensus (row) vs all consensuses above
# Thus, this is particularly programmed to represent a "consensus similarity gorup" by the consensus (row)
#     with the highest p-value (i.e. collapse and sort by p-value)
# Takes: a list of Consensus Motifs
# Returns: A dataframe with motifs as rows, and second col as max levenshteinSim value with rest of list
helper_getMaxLevenSimCol <- function(vectorOfMotifs){
  # Convert vector to dataframe
  vectorDF <- as.data.frame(vectorOfMotifs)
  # Make the levenSim of row 1 equal to 0, since it has no rows above to compare. This simplifies later code
  vectorDF[1,2]=0.0
  # Record how many rows/motifs to process, used for for loop
  vectorLen <- length(vectorOfMotifs)
  # For loop to go through each row, compare the row's consensus to all previous consesuses, and store as a column (col2)
  for(i in 2:vectorLen){
    vectorDF[i,2]=max(sapply(as.list(vectorOfMotifs[1:(i-1)]), function(x) RecordLinkage::levenshteinSim(x, vectorOfMotifs[i])))
  }
  # Rename columns of this output dataframe
  colnames(vectorDF) <- c("Consensus","LevenSim")
  # Return
  return(vectorDF)
}

