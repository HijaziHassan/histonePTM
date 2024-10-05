
#' Split spectrum_title column into scan and file columns
#'
#' @param df dataframe
#' @param col The column to be fragmented.
#' @param scan_prefix The pattern the preceeds the scan number (e.g. 'Index:' or "Scan:"). If there is a space before it make sure to include it.
#' @param file_prefix The pattern the preceeds the file name (e.g. 'raw_file:). If there is a space before it make sure to include it.
#' @param sep A separator as ";" or ":" (`optional`).
#' @param ... to add more arguments that could be useful as this function is a wrapper around `separate_wider_regex` function.
#'
#' @importFrom tidyr separate_wider_regex
#' @importFrom cli cli_alert_info
#' @return two columns: `id_scan` and `id_filename`
#' @export
#'

misc_parseSpectrumtitle <- function(df, col,  scan_prefix, file_prefix, sep ="", ...){



  if(!{{col}} %in% colnames(df)) {

    cli::cli_alert_info("{col} column has either been already parsed or doesn't exist")

    }else{


  regex_any=  '.*?'
  scan_regex = '\\d+'
  file_regex = '.+'

  patterns = c(regex_any, scan_prefix, id_scan = scan_regex, regex_any, file_prefix, id_filename = file_regex, paste0(sep, regex_any))

df <- df |>
  tidyr::separate_wider_regex(cols = {{col}},
                              patterns = patterns,
                              cols_remove = TRUE,
                              too_few = "align_start",
                              ...)

}

return(df)



}


