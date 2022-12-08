#' Plot HOMER TF enrichment results
#'
#' # Strategy:
#' ### Find motifs to extract
#' (1) Filter each element table by q-value. Should basically chop off the bottom portion of list, making the rest of the script less computationally cumbersome
#' (2) Collapse each element table by by Levenshtein Similarity
#' (3) Filter each element table (in my case: cell type-specific results file) to top X rows
#' (4) Extract consensus columns from each element table and store as a variable]
#'
#' ### Make bigTable of all TFs to pull from so a single TF can have data from different input Files (e.g. across CTs)
#' (5) Create bigTable of all q-value TF results concatenated together
#' (6) Filter by consensus list (4)and make a gg-plot appropriate table and Plot
#' (7-8) Order factors and plot

#'
#' @param dir string. Input directory containing HOMER findMotifsGenome.pl output files in format: *knownResults.txt
#' @param show int. Number of rows to show per input file, ranked by p-value.
#' @param qThreshold int. Value for thresholding HOMER enrichment results by q-value.
#' @param levenSimThreshold float. Value for thresholding TFs. For groups of TFs with similar consensus sequences, the TF with the lowest p-value by HOMER will be retained.
#'
#' @author Tim Scott
#'
#' @return ggplot object
#'
#' @export
plot_HOMERTFs <- function(dir="/directory/of/results/", show=3, qThreshold = 0.05, levenSimThreshold=1) {
    ### Import files:
    # Find names and paths of files matching standard "*knownResults.txt" pattern from results directory
    loFilenames = list.files(path = dir, pattern = "*knownResults.txt$")
    loFilepaths = paste0(dir, loFilenames)
    # Use list of path names to load in files
    loFiles = lapply(loFilepaths, readr::read_tsv, col_names=c("MotifName","Consensus","p","logp","q","numTarget","percentTarget","numBackground","percentBackground"), skip = 1)

    ### Clean up imported files:
    # Use lapply to rename each table in the list based on the filename as in the directory
    # Use the first list of filenames we used to load in to grab the filenames
    loFilenames_short <- lapply(loFilenames, function(x) gsub("*.knownResults.txt", "", x))
    names(loFiles) <- loFilenames_short
    loFiles_pFilt <- lapply(loFiles, function(x) dplyr::filter(x, .data$p<0.05))
    loFiles_Enrich <- lapply(loFiles_pFilt, function(x) {
    dplyr::mutate(x, percentTargetNum = (as.numeric(gsub("%", "", .data$percentTarget, fixed = TRUE))/100),
                      percentBackgroundNum = (as.numeric(gsub("%", "", .data$percentBackground, fixed = TRUE))/100),
                      percentFold = (.data$percentTargetNum/.data$percentBackgroundNum))})


    ### Process the input element tables to determine a list of TF consesunses we want to display
    # (1) Filter each input list element by q-value
    loFiles_Enrich_qValFilt <- lapply(loFiles_Enrich, function(x) dplyr::filter(x, q <= qThreshold))
    # (2) Collapse each input list element by Levenshtein Similarity threshold
    loFiles_Enrich_qValFilt_LevenSimFilt <- lapply(loFiles_Enrich_qValFilt, helper_collapseTableByLevenSim,
                                                   levenSimThresholdVal=levenSimThreshold)
    # (3) Extract the top X rows of TFs from each input list element (e.g. per each Celltype)
    loFiles_Enrich_qValFilt_LevenSimFilt_topX <- lapply(loFiles_Enrich_qValFilt_LevenSimFilt, '[', c(1:show),)
    # (4) Get df of processed/filtered concensuses
    dfOfFilteredConcensuses <- data.table::rbindlist(lapply(loFiles_Enrich_qValFilt_LevenSimFilt_topX, function(x) x[,2]))


    ### Make a bigTable of all results to pull from
    # (5) Cat and add new cols: Clean % numerical values and motif names (filtered "MotifName"before first open paren)
    bigTable <- helper_makeBigTableFromListOfStandardTables(loFiles_Enrich_qValFilt)
    # (6) Pull rows from bigTable (5) as defined in (4)
    tableToPlot <- dplyr::filter(bigTable, .data$Consensus %in% dfOfFilteredConcensuses[[1]])
    # (7) Order factors for ggplot
    # Change to ID (e.g. Celltype) factor to maintain order in plotting. Please customize.
    tableToPlot$ID <- factor(tableToPlot$ID,
                             levels=loFilenames_short)
    # This command reorders X-axis motifs nicely.
    tableToPlot$Motif <- factor(tableToPlot$Motif, levels=rev(unique(tableToPlot$Motif[(order(tableToPlot$ID))])))


    ### Plot with GGPlot
    # Find the maximum P value to scale the color with
    maxPValue <- max(abs(tableToPlot$logp))
    minPValue <- min(abs(tableToPlot$logp))
    # (8) Set colors and plot
    col.pal <- RColorBrewer::brewer.pal(6,"YlOrRd")

    # (9) Call ggplot function
    hm <- ggplot2::ggplot(tableToPlot, ggplot2::aes(x=.data$Motif, y=.data$ID)) +
      # ggplot2::geom_tile(aes(fill = -(logp)),colour="white",size=0.5) + # size=.5 adds a little white spacing between tiles
      ggplot2::geom_point(ggplot2::aes(size = .data$percentFold, color = abs(.data$logp))) +
      ggplot2::scale_colour_gradient2(limits=c(minPValue, maxPValue), low="yellow", mid="orange", high="red") +
      ggplot2::ggtitle("HOMER Motif") +
      ggplot2::scale_size_continuous(range = c(3,8)) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
            axis.text.x = ggplot2::element_text(angle=-90, hjust=0, vjust=0),
            axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::scale_x_discrete(position = "bottom")  + # this moves the x-axis to the top, from the bottom
      ggplot2::scale_y_discrete(position = "left") +
      ggplot2::ylab(NULL) +
      ggplot2::theme(legend.position="right") # needed horizontal space personally, so this moves legend to bottom
    hm # return heatmap
} # End of function definition


