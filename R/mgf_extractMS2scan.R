


#function to read mgf file
.mgf_to_sp <- function(mgf_file){
  cli::cli_inform(message = "Converting {mgf_file} file into Spectra object")

  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}

#extract scan to mgf
.scan_to_df <- function(scan, spec, mgf_file){

  spec <- spec[spec$acquisitionNum == scan]

  todf <- data.frame(
    scan = scan,
    mgf= basename(mgf_file),
    unlist(peaksData(spec))
  )


}



#' Extract mz and intensity values of an MS2 scan from an mgf file
#' @description
#' extract user-defined scans as \code{mgf} or as \code{.csv} file.
#'
#' @param mgf_file .mgf file.
#' @param scan scan number(s) containing of MS/MS events of interest.
#' @param save_file \code{logical}. Save scans' mz and intensity values as \code{.csv} file.
#' @param export_mgf \code{logical}. export subset mgf file only containing the selected MS/MS scans.
#'
#' @return \code{dataframe} containing 4 columns: _mz_, _intensity_, _scan_, and _file_.
#'
#' @import Spectra
#' @import MsBackendMgf
#' @import purrr
#' @import BiocParallel
#' @import cli
#' @import stringr
#'
#'
#'@export
mgf_extractMS2scan <- function(mgf_file, scan,  save_file = FALSE, export_mgf = FALSE){
  if(rlang::is_missing(mgf_file)) cli::cli_abort(c("Error:",
                   "i" = 'argument mgf_file is missing with no default'))

  file_name = basename(stringr::str_remove(string = mgf_file, pattern = ".mgf"))


        sps <- .mgf_to_sp(mgf_file)

    sub_sps <- sps[sps$acquisitionNum %in% c(scan)]


    if(isTRUE(length(sub_sps$acquisitionNum) == 0L)){
      cli::cli_abort(c("Must provide valid scan number(s).",
                       "i"= " The following {scan} do(es)n't exist in {mgf_file} file." ,
                       "x"= "You can quickly check by opening {mgf_file} in any text editor."))
    }



     if (isTRUE(export_mgf)) {
    file_submgf = paste0("selectedscans_", file_name, ".mgf")
    map <-
      c(spectrumName = "TITLE", Spectra::spectraVariableMapping(MsBackendMgf::MsBackendMgf()))
    Spectra::export(sub_sps,
           backend = MsBackendMgf::MsBackendMgf(),
           file = file_submgf,
           mapping = map)
    cli::cli_alert_success("{file_submgf} file is saved sucessfully.")
  }



  sub_sps_df <- map(.x = scan, ~ .scan_to_df(
    scan = .x,
    spec = sub_sps,
    mgf_file = mgf_file
  )) |>
    dplyr::bind_rows() |>
    as_tibble()


  return(sub_sps_df)

  if (isTRUE(save_file)) {
    file_csv = paste0("selectedscans_", file_name, ".csv")
    write.csv(x =  sub_sps_df,
              file = file_csv,
              row.names = FALSE)
    cli::cli_alert_success("{file_csv} file is saved sucessfully.")
  }



}


