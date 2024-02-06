

#function to read mgf file
.mgf_to_sp <- function(mgf_file){
  Spectra::Spectra(object = {{mgf_file}}, source = MsBackendMgf::MsBackendMgf())

}



#' Find MS/MS spectra that contain a diagnostic ion.
#'
#' @param mgf_file mgf file to search.
#' @param diag_ion  The mz of the diagnsotic ion
#' @param tol A mass tolerance to respect during the search
#'
#' @return A \code{tibble} with 6 columns.
#' @import dplyr
#' @import Spectra
#' @import MsBackendMgf
#'
diagIon_find <- function(mgf_file, diag_ion, tol = 0.01){

  sps <- mgf_to_sp(mgf_file)

  list = list()

  for (i in 1:length(acquisitionNum(sps))){
    isIonhere <-  dplyr::near(diag_ion, mz(sps)[[i]], tol = tol) #logical vector
    if(isTRUE(any(isIonhere, na.rm = TRUE))) #find at least one is TRUE

      indx =   which(isIonhere) #find its index

    #report data
    list[[i]] <- c(file      = mgf_file,
                   diag_ion  = diag_ion,
                   int_diag  = Spectra::intensity(sps)[[i]][indx],
                   scan      = Spectra::acquisitionNum(sps)[[i]],
                   prec_mass = Spectra::precursorMz(sps)[[i]],
                   prec_z    = Spectra::precursorCharge(sps)[[i]])
  }

  bind_rows(list)
}


