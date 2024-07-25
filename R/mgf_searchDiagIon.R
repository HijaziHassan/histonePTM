
#' Find MS/MS spectra that contain a diagnostic ion.
#' @description
#' Report user-defined diagnsotic ion(s) per MS/MS spectrum from an \code{mgf} file.
#' @param mgf_file mgf file to search.
#' @param diag_ion The mz of the diagnsotic ion
#' @param tol A mass tolerance to respect during the search
#' @param save_file \code{logical} Save the results as csv file.
#'
#' @return A \code{tibble} with 6 columns including the diagnostic ion \code{m/z} and its intensity relative to the base peak.
#' @import dplyr
#' @import Spectra
#' @importFrom MsBackendMgf MsBackendMgf
#' @import Spectra
#' @import cli
#' @import rlang
#' @import BiocManager
#'
#' @export
mgf_searchDiagIon <- function(mgf_file, diag_ion, tol = 0.002, save_file = FALSE){






  libraries <- c("BiocParallel", "Spectra", "MsBackendMgf")

  if (!requireNamespace("BiocManager")) install.packages("BiocManager")

  #Checking if the package belongs to CRAN or to Bioconductor and installing them accordingly.

  for(lib in libraries){
    if(!lib %in% installed.packages()){
      if(lib %in% available.packages()[,1]){
        install.packages(lib,dependencies=TRUE)
      }else {(BiocManager::install(lib))
      }}
  }

  #Loading the libraries
  sapply(libraries,requireNamespace,character=TRUE)


if(rlang::is_missing(mgf_file)){
  cli::cli_abort("No mgf file name is identified.")
}
if(is_missing(diag_ion)){
  cli::cli_abort("No diagnostic ion is identified.")
}



  file_name = tools::file_path_sans_ext(basename(mgf_file))

  #read mgf file(s) as Spectra object(s)
  sps <- .mgf_to_sp(mgf_file)


  final_list = vector(mode = 'list')


#cli::cli_inform(message = "Extracting diagnostic ions from {mgf_file}.\n\n")

cli::cli_progress_bar(type = "iterator", name = paste0("Extracting diagnostic ions from ", mgf_file))

  for (i in 1:length(Spectra::acquisitionNum(sps))) {

    temp_all_diag <- list()


    for (ion in seq_along(diag_ion)) {

      isIonhere <- dplyr::near(diag_ion[ion], Spectra::mz(sps)[[i]], tol = tol) #logical vector

      #find at least one is TRUE
      if (any(isIonhere, na.rm = TRUE)) {
        indices =   which(isIonhere) #find its/their index/indices

        temp_one_diag <- list() # to account for two values too close to each other.

        for (indx in seq_along(indices)) {
          #report data
          temp_one_diag[[indx]] <- c(
            file         = file_name,
            diag_ion     = round(Spectra::mz(sps)[[i]][indices[indx]], 4),
            diag_relint  = round(Spectra::intensity(sps)[[i]][indices[indx]] / max(Spectra::intensity(sps)[[i]]),3),
            #scan         = Spectra::acquisitionNum(sps)[[i]],
            scan         = .found_replace_scan(sps, i),
            prec_mz      = round(Spectra::precursorMz(sps)[[i]], 4),
            prec_z       = Spectra::precursorCharge(sps)[[i]],
            rt           = round(Spectra::rtime(sps)[[i]] / 60, 2) #sec to min
          )

          cli::cli_progress_update()

        }


        # Update the combined list
        temp_all_diag <- append(temp_all_diag, temp_one_diag)


      }


    }


    final_list[[i]] <- dplyr::bind_rows(temp_all_diag)
  }

  final_df <- dplyr::bind_rows(final_list)

  if(nrow( final_df) == 0L){cli::cli_inform("No diagnostic ion was found.")}else{

    final_df <- final_df |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("diag_ion", "diag_relint", "prec_mz", "rt")), as.double),
      dplyr::across(dplyr::all_of(c("scan", "prec_z")), as.integer)
    )



    return(final_df)


    if(save_file){
      file_csv = paste0("diagIons_", file_name, ".csv")
      write.csv(x = final_df, file = file_csv, row.names = FALSE)
      wd = getwd()
      cli::cli_alert_success("{file_csv} file is saved sucessfully into {wd}.")
    }

}
cli::cli_progress_done(result = "done")







}

#' @examples
#' # example code
#'



#' @noRd
.mgf_to_sp <- function(mgf_file){

  BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
  mgftospec <- Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())
  mgf_file_basename = basename(mgf_file)
cli::cli_alert_success( "The file {mgf_file_basename} is converted into a 'Spectra' object.")
return(mgftospec)
}

#' @noRd

.found_replace_scan= function(spec_obj, scan_i){

  ifelse(!is.na(spec_obj$acquisitionNum[scan_i]),
         spec_obj$acquisitionNum[scan_i],
         {

           title = spec_obj[['TITLE']][scan_i]

           pattern <- "(scan=|_?scan:|Scan|Index:)\\s*(\\d+)"

           # Search for the pattern in the title
           match <- stringr::str_extract(title, pattern = pattern, group = 2)
           print(match)

           #if scan is found in the title, replace NA with this scan to extract mz and intensity values
           if(!is.na(match)){
             scan_number <- as.integer(match)
            spec_obj$acquisitionNum[scan_i] = scan_number
           }else{NA}
         }
  )
}

