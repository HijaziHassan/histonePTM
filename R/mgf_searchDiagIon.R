


.mgf_to_sp <- function(mgf_file){
  cli::cli_inform(message = "Converting {mgf_file} file into Spectra object")
  BiocParallel::register(BiocParallel::SerialParam())
  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}



#' Find MS/MS spectra that contain a diagnostic ion.
#' @description
#' Report user-defined diagnsotic ion(s) per MS/MS spectrum from an \code{mgf} file.
#' @param mgf_file mgf file to search.
#' @param diag_ion The mz of the diagnsotic ion
#' @param tol A mass tolerance to respect during the search
#' @param save_file(logical) Save the results as csv file.
#'
#' @return A \code{tibble} with 6 columns including the diagnostic ion \code{m/z} and its intensity relative to the base peak.
#' @import dplyr
#' @import Spectra
#' @import BiocParallel
#' @import MsBackendMgf
#' @import cli
#' @export
mgf_searchDiagIon <- function(mgf_file, diag_ion, tol = 0.002, save_file = FALSE){

  file_name = basename(stringr::str_remove(string = mgf_file, pattern = ".mgf"))

  sps <- .mgf_to_sp(mgf_file)


  final_list = list()


cli::cli_inform(message = "Extracting diagnostic ions from {mgf_file}.")

cli::cli_progress_bar()

  for (i in 1:length(acquisitionNum(sps))) {
    temp_list_combined <- list()


    for (ion in seq_along(diag_ion)) {
      isIonhere <-
        dplyr::near(diag_ion[ion], mz(sps)[[i]], tol = tol) #logical vector

      #find at least one is TRUE
      if (isTRUE(any(isIonhere, na.rm = TRUE))) {
        indices =   which(isIonhere) #find its/their index/indices

        temp_vec_list <- list() # to account for two values too close to each other.

        for (indx in seq_along(indices)) {
          #report data
          temp_vec_list[[indx]] <- c(
            file         = file_name,
            diag_ion     = round(Spectra::mz(sps)[[i]][indices[indx]], 4),
            diag_relint  = round(Spectra::intensity(sps)[[i]][indices[indx]] / max(Spectra::intensity(sps)[[i]]),3),
            scan         = Spectra::acquisitionNum(sps)[[i]],
            prec_mass    = round(Spectra::precursorMz(sps)[[i]], 4),
            prec_z       = Spectra::precursorCharge(sps)[[i]],
            rt           = round(Spectra::rtime(sps)[[i]] / 60, 2) #sec to min
          )

          cli::cli_progress_update()

        }

        # Update the combined list
        temp_list_combined <- c(temp_list_combined, temp_vec_list)


      }


    }


    final_list[[i]] <- dplyr::bind_rows(temp_list_combined)
  }

  final_df <- dplyr::bind_rows(final_list) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("diag_ion", "diag_relint", "prec_mass", "rt")), as.double),
      dplyr::across(dplyr::all_of(c("scan", "prec_z")), as.integer)
    )

#cli::cli_progress_done()

  return(final_df)

  if(isTRUE(save_file)){
    file_csv = paste0("diagIons_", file_name, ".csv")
    write.csv(x = df, file = file_csv, row.names = FALSE)
    cli::cli_alert_success("{file_name} file is saved sucessfully.")
  }

}

#' @examples
#' # example code
#'


