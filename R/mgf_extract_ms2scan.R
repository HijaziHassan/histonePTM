

#function to read mgf file
.mgf_to_sp <- function(mgf_file){
  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}


#function to extract scan from Spectra object created out of mgf
.extract_scan <- function(sps, scan, mgf_file){

  peaks_df <- sps[sps$acquisitionNum %in% c(scan)] |>
    peaksData() |>
    unlist()|>
    as.data.frame()

  peaks_df$scan <- scan
  peaks_df$file <- mgf_file

  return(peaks_df)

}


#' Extract mz & intensity of an MS2 scan from an mgf file
#'
#' @param mgf_file .mgf file.
#' @param scan scan number containing ms2 of interest
#'
#' @return \code{dataframe} containing 4 columns: _mz_, _intensity_, _scan_, and _file_.
#'
#' @import Spectra
#' @import MsBackendMgf
#' @import purrr
#'
#'
#'@export
mgf_extract_ms2scan <- function(mgf_file, scan){
  sps <- .mgf_to_sp(mgf_file)

  df <- purrr::map_df(scan, .f = ~.extract_scan(sps = sps, scan= .x, mgf_file = mgf_file))

  return(df)
}



