


.mgf_to_sp <- function(mgf_file){

  BiocParallel::register(BiocParallel::SerialParam())

  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}



#' Find MS/MS spectra that contain a diagnostic ion.
#' @description
#' Report user-defined diagnsotic ion(s) per MS/MS spectrum.
#' @param mgf_file mgf file to search.
#' @param diag_ion  The mz of the diagnsotic ion
#' @param tol A mass tolerance to respect during the search
#'
#' @return A \code{tibble} with 6 columns diagnostic ion and its intensity relative to the base peak.
#' @import dplyr
#' @import Spectra
#' @import BiocParallel
#' @import MsBackendMgf
#' @import cli
#' @export
mgf_searchDiagIon <- function(mgf_file, diag_ion, tol = 0.002){

  cli::cli_inform(message = "Converting {mgf_file} file into Spectra object")


  sps <- .mgf_to_sp(mgf_file)


  final_list = list()


  cli::cli_inform(message = "Extracting diagnostic ions from {mgf_file}.")


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
            file         = mgf_file,
            diag_ion     = round(Spectra::mz(sps)[[i]][indices[indx]], 4),
            diag_relint  = round(Spectra::intensity(sps)[[i]][indices[indx]] / max(Spectra::intensity(sps)[[i]]),3),
            scan         = Spectra::acquisitionNum(sps)[[i]],
            prec_mass    = round(Spectra::precursorMz(sps)[[i]], 4),
            prec_z       = Spectra::precursorCharge(sps)[[i]],
            rt           = round(Spectra::rtime(sps)[[i]] / 60, 2) #sec to min
          )



        }

        # Update the combined list
        temp_list_combined <- c(temp_list_combined, temp_vec_list)


      }


    }


    final_list[[i]] <- dplyr::bind_rows(temp_list_combined)
  }

  final_df <- dplyr::bind_rows(final_list)



  return(final_df)


}

#' @examples
#' # example code
#'


