




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
#' @importFrom cli cli_abort cli_inform cli_alert_success
#' @importFrom stringr str_remove str_extract
#' @importFrom utils write.csv
#'
#'
#'@export
mgf_extractMS2scan <- function(mgf_file, scan,  save_file = FALSE, export_mgf = FALSE){
  if (!requireNamespace("BiocManager")) install.packages("BiocManager")

  if (!requireNamespace("Spectra", quietly = TRUE)) {
    cli::cli_abort("Package 'Spectra' is required but not installed. Please install it using BiocManager::install('Spectra').")
  }

  if (!requireNamespace("MsBackendMgf", quietly = TRUE)) {
    cli::cli_abort("Package 'MsBackendMgf' is required but not installed. Please install it using BiocManager::install('MsBackendMgf').")
  }

  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    cli::cli_abort("Package 'BiocParallel' is required but not installed. Please install it using BiocManager::install('BiocParallel').")
  }


  if(rlang::is_missing(mgf_file)) cli::cli_abort(c("Error:",
                   "i" = 'argument mgf_file is missing with no default'))

  file_name = basename(stringr::str_remove(string = mgf_file, pattern = ".mgf"))


        sps <- .mgf_to_sp(mgf_file)

    sub_sps <- sps[sps$acquisitionNum %in% c(scan)]

    if(length(sub_sps$acquisitionNum) == 0L & all(is.na(sps$acquisitionNum))){

      pattern <- "(scan=|_?scan:|Scan|Index:)\\s*(\\d+)"
      for(s in scan){
      for(ti in 1:length(sps)){


        title = sps$TITLE[ti]

        # Search for the pattern in the title
        match <- stringr::str_extract(title, pattern = pattern, group = 2)

        #if scan is found in the title, replace NA with this scan to extract mz and intensity values
        if(match == s){
          scan_number <- as.integer(match)
          sps$acquisitionNum[ti] = scan_number
          cat("Scan number found:", scan_number, "\n")
          break

        }else{next}

        }
      }

      sub_sps <- sps[sps$acquisitionNum %in% c(scan)]

      }else{

      cli::cli_abort(c("Must provide valid scan number(s).",
                       "i"= " The following {scan} do(es)n't exist or it is neither found
                       as a separate line or in the TITLE line in {mgf_file} file." ,
                       "x"= "You can quickly check by opening {mgf_file} in any text editor."))
    }



     if(export_mgf){
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
    tidyr::as_tibble()


  return(sub_sps_df)

  if (save_file) {
    file_csv = paste0("selectedscans_", file_name, ".csv")
    write.csv(x =  sub_sps_df,
              file = file_csv,
              row.names = FALSE)
    cli::cli_alert_success("{file_csv} file is saved sucessfully.")
  }



}



#' @noRd

#function to read mgf file
.mgf_to_sp <- function(mgf_file){
  cli::cli_inform(message = "Converting {mgf_file} file into Spectra object")
  BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}



#' @noRd
#extract scan to mgf
.scan_to_df <- function(scan, spec, mgf_file){

  spec <- spec[which(spec$acquisitionNum == scan)]

  todf <- data.frame(

    mgf= basename(mgf_file),
    scan = as.integer(scan),
    unlist(peaksData(spec))

  )


}

