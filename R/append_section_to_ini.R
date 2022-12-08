#' Append section to ini file
#'
#' Takes a new section in ini format and adds to existing ini.
#'
#' The new_section must be a named list of the section list. See examples.
#'
#'
#' @param ini_file file location of config.ini file
#' @param new_section named list of the section list
#' @import ini
#' @author Tyler Hansen
#' @examples
#' #list of key-value pairs
#' CHRACC <- list(dir='/chrAcc_peaks/',
#'                peaks='/chrAcc_peaks/GM12878_genrich.narrowPeak')
#'
#' #list of section, resulting in list of list.
#' new_section <- list(CHRACC=CHRACC)
#'
#' #write ini
#' ini_file <- system.file("extdata", "config.ini")
#' append_section_to_ini(ini_file, new_section)
#'
#' @return No return value. Edits and overwrites input config.ini file.
#'
#' @export
append_section_to_ini <-
function(ini_file, new_section) {
    if (is.list(new_section)) {
        ini <- ini::read.ini(ini_file)
        ini <- c(ini,new_section)
        ini::write.ini(x=ini, filepath=ini_file)
    } else {
        stop('new_section is not a list')
    }
}
