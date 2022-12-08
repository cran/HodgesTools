#' CpG Analysis
#'
#' Compute observed/expected CpG ratio and GC% for regions of interest
#'
#' The function reads in a nucOutput_gc and CpGcount text file
#' The function uses the nucOutput_gc and CpGcount file to calculates observed/expected ratio and GC%.
#' The function allows the option to plot the distribution of these values in ggplot2
#'
#' @param list "boolean of whether input is a list of groups. Default = FALSE."
#' @param cpg_file "file names or list of files names for your CpGcount.txt files. This is defined in cpg_analysis.sh"
#' @param nuc_file "file names or list of files names for your nucOutput_gc.txt files. This is defined in cpg_analysis.sh"
#' @param count "numeric value for the number of files included in your list
#' @param palette "if choosing to plot, the RColorBrewer palette you would like to be applied to your plot"
#' @param plot "one of three choices depending on what output you would like: 'none' for no plot, 'ratio' for observed/expected ratios, 'gc_percent' for GC%"
#'
#' @import tools
#'
#' @author Lindsey Guerin
#'
#' @export
#'
#' @return ggplot object or tibble if plot="none"
#'
#' @examples
#' #load external data
#'
#' gain_6hr_CpG <- system.file(package = "HodgesTools", "extdata",
#' "cov5root_6hr_gain.CpGcount.txt")
#' gain_12hr_CpG <- system.file(package = "HodgesTools", "extdata",
#'  "cov5root_12hr_gain.CpGcount.txt")
#' gain_6hr_nuc <- system.file(package = "HodgesTools", "extdata",
#'  "cov5root_6hr_gain.nucOutput_gc.txt")
#' gain_12hr_nuc <- system.file(package = "HodgesTools", "extdata",
#' "cov5root_12hr_gain.nucOutput_gc.txt")
#'
#' #Make a density plot of GC% values for a list of two region of interest files
#' cpg_analysis(list = TRUE, count = 2, cpg_file = list(gain_6hr_CpG, gain_12hr_CpG),
#' nuc_file= list(gain_6hr_nuc, gain_12hr_nuc), palette = "Set3", plot ="gc_percent")
#'
#' #Make a density plot of observed/expected values for a single set of regions of interest
#' cpg_analysis(list = FALSE, cpg_file = gain_6hr_CpG,
#' nuc_file = gain_6hr_nuc, palette = "Set3", plot ="ratio")
cpg_analysis <- function(list = FALSE, count, cpg_file, nuc_file, palette = "Set3", plot = "none"){
  if (list == FALSE){
    #read in files and assign colnames
    cpgcount <- readr::read_tsv(cpg_file, col_names = c("CpGs", "chr", "start", "end"), col_types = "ncnn")
    nucOutput <- readr::read_tsv(nuc_file, col_names = c("chr", "start", "end", "gc%", "numC", "numG", "length"), col_types = "cnnnnnn")
    #join obeserved and expected data and calculate values
    join <- dplyr::inner_join(cpgcount, nucOutput)
    join <- dplyr::mutate(join, ratio= (.data$CpGs*.data$length)/(.data$numC*.data$numG))

    if (plot == "ratio"){
      ratio_plot <- ggplot2::ggplot(join, ggplot2::aes(x=.data$ratio, alpha= 0.5))+
        ggplot2::geom_density()+
        ggplot2::scale_color_brewer(palette= palette)+
        ggplot2::scale_fill_brewer(palette= palette)+
        ggplot2::theme_minimal()+
        ggplot2::xlim(0,1)+
        ggplot2::xlab("Obs/Exp CpG ratio")+
        ggplot2::ylab("Proportion")+
        ggplot2::ggtitle("Obs/Exp CpG Ratio by Group")
      ratio_plot
    }else if (plot == "gc_percent"){
      gccontent_plot <- ggplot2::ggplot(join, ggplot2::aes(x=.data$`gc%`, alpha= 0.5))+
        ggplot2::geom_density()+
        ggplot2::scale_color_brewer(palette= palette)+
        ggplot2::scale_fill_brewer(palette= palette)+
        ggplot2::theme_minimal()+
        ggplot2::xlab("GC%")+
        ggplot2::ylab("Proportion")+
        ggplot2::ggtitle("GC% by Group")
      gccontent_plot

    }else if (plot =="none"){
      message("Not plotting results")
      join
    }
  } else if (list == TRUE){

    for(i in 1:count){
      regions <-sub("\\.CpGcount.txt", "", basename(cpg_file[[i]]))
      cpgcount <- readr::read_tsv(cpg_file[[i]], col_names = c("CpGs", "chr", "start", "end"), col_types = "ncnn")
      nucOutput <- readr::read_tsv(nuc_file[[i]], col_names = c("chr", "start", "end", "gc%", "numC", "numG", "length"), col_types = "cnnnnnn")

      #join observed and expected data and calculate values
      singleGroup_join <- dplyr::inner_join(cpgcount, nucOutput)
      singleGroup_join <- dplyr::mutate(singleGroup_join, ratio= (.data$CpGs*.data$length)/(.data$numC*.data$numG))
      singleGroup_join <- dplyr::mutate(singleGroup_join, group = regions)
      singleGroup_join
    }
    singleGroup_join$group<- as.factor(singleGroup_join$group)

    if (plot == "ratio"){
      ratio_plot <- ggplot2::ggplot(singleGroup_join, ggplot2::aes(x=.data$ratio, fill=.data$group, color = .data$group, alpha= 0.5))+
        ggplot2::geom_density()+
        ggplot2::facet_wrap(singleGroup_join$group)+
        ggplot2::scale_color_brewer(palette= palette)+
        ggplot2::scale_fill_brewer(palette= palette)+
        ggplot2::theme_minimal()+
        ggplot2::xlim(0,1)+
        ggplot2::xlab("Obs/Exp CpG ratio")+
        ggplot2::ylab("Proportion")+
        ggplot2::ggtitle("Obs/Exp CpG Ratio by Group")
      ratio_plot
    }else if (plot == "gc_percent"){
      gccontent_plot <- ggplot2::ggplot(singleGroup_join, ggplot2::aes(x=.data$`gc%`, fill=.data$group, color = .data$group, alpha= 0.5))+
        ggplot2::geom_density()+
        ggplot2::facet_wrap(singleGroup_join$group)+
        ggplot2::scale_color_brewer(palette= palette)+
        ggplot2::scale_fill_brewer(palette= palette)+
        ggplot2::theme_minimal()+
        ggplot2::xlab("GC%")+
        ggplot2::ylab("Proportion")+
        ggplot2::ggtitle("GC% by Group")
      gccontent_plot

    }else if (plot =="none"){
      message("Not plotting results")
      singleGroup_join
    }
  }
}
