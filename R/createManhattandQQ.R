#' Creating a Manhattan Plot and QQ plot
#'
#' Creates a Manhattan plot and QQ plot using GWAS results output from PLINK
#'
#' This function reads in a GWAS result file output from plink2 listing the coordinates, ids, and associated p-values for SNPs under study
#' This function also has the option of reading in a "highlights" file listing the IDs of SNPs to annotate/color on the Manhattan plot
#'
#' @param gwas_results output file listing SNP-trait association values for GWAS run using PLINK
#' @param highlights_file a text file with a 'snp' column listing the SNPs to annotate/color on the Manhattan plot
#' @param suggestive_line where to draw a "suggestive" line; default -log10(1e-5).
#' @param genomewide_line where to draw a "genome-wide significant" line; default -log10(5e-8)
#' @param annotate_Pval if set, SNPs below this p-value will be annotated on the plot; default is 0.05
#' @param y_lim set the y-axis limits; default is c(0,8)
#' @param set_color_vector a character vector listing colors in palette of interest (you must create this chr object before calling the createManhattanandQQ function and assign it to set_color_vector)
#'
#' @return a Manhattan plot of SNP-trait associations and QQ plot
#'
#' @import qqman
#' @import dplyr
#'
#' @author Verda Agan
#' @examples
#' #' #load external data.
#' gwas_results <- system.file(package = "HodgesTools", "extdata",
#' "createManhattandQQ_example_sum_stats.txt")
#' snps_to_annotate <- system.file(package = "HodgesTools", "extdata",
#' "createManhattandQQ_example_highlights_file.txt")
#'
#' #Make a Manhattan plot that highlights a select list of SNPs subset from GWAS results
#' createManhattandQQ(gwas_results, highlights_file=snps_to_annotate,
#' suggestive_line = -log10(0.001), set_color_vector = c("gray10", "gray60"),
#' genomewide_line = -log10(5e-8), annotate_Pval = 0.001, y_lim =c(0,8))
#'
#' #Make a Manhattan plot that doesn't highlight a select list of SNPs subset from GWAS results
#' createManhattandQQ(gwas_results, suggestive_line = -log10(0.001),
#' set_color_vector = c("gray10", "gray60"), genomewide_line = -log10(5e-8),
#' annotate_Pval = 0.001, y_lim =c(0,8))
#' @export

#col_vector <- brewer.pal(n = 8, name = "Paired")
createManhattandQQ <- function(gwas_results, highlights_file=NULL, suggestive_line = -log10(0.05), set_color_vector = c("gray10", "gray60"), genomewide_line = -log10(5e-8), annotate_Pval = 0.05, y_lim = c(0,8)) {
###set variables
message(paste0("suggestive line is set to: ", suggestive_line))
message(paste0("genomewide line is set to: ", genomewide_line))
message(paste0("annotatePval is set to: ", annotate_Pval))

gwas_results <- readr::read_tsv(file = gwas_results, col_names = TRUE)
gwas_results <- gwas_results %>% dplyr::select("#CHROM", "POS", "ID", "P") %>% magrittr::set_colnames(c("chrom","start","snp","pval"))

#this if statement will parse the function depending on whether a highlights_file is specified
if (!is.null(highlights_file)) { #TRUE
  #read in file listing SNPs of interest to be annotated on Manhattan plot
  snps_to_highlight <- readr::read_table(file = highlights_file, col_names = TRUE)
  snps_to_highlight <- snps_to_highlight %>% dplyr::select("snp")
  gwas_results_for_subset_highlight_snps <- gwas_results %>% dplyr::filter(.data$snp %in% snps_to_highlight$snp)
  #make manhattan plot
  manhattan_plot <- manhattan(gwas_results, chr="chrom", bp="start", p="pval", snp="snp",
                              main = "Manhattan plot: Logistic", ylim=y_lim, logp = TRUE, annotatePval = annotate_Pval, col = set_color_vector,
                              highlight = gwas_results_for_subset_highlight_snps$snp, annotateTop = FALSE,
                              suggestiveline = suggestive_line, genomewideline = genomewide_line)
  } else { #FALSE
    #if highlight parameter is FALSE because we don't want to highlight select SNPs on the Manhattan
    #read in gwas results and assign column names
    #make manhattan plot
    manhattan_plot <- manhattan(gwas_results,chr="chrom",bp="start",p="pval",snp="snp",
                                main = "Manhattan plot: Logistic", ylim=y_lim, logp = TRUE, annotatePval = annotate_Pval, col = set_color_vector,
                                annotateTop = FALSE,
                                suggestiveline = suggestive_line, genomewideline = genomewide_line)
  }

#make qqplot
qq_plot <- qq(gwas_results$pval, main = "Q-Q plot of GWAS p-values")

#return objects
return(list(MP=manhattan_plot,QQ=qq_plot))
}
