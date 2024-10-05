
#' Find MS/MS spectra that contain a diagnostic ion.
#' @description
#' Report user-defined diagnsotic ion(s) per MS/MS spectrum from an \code{mgf} file.
#' @param mgf_file mgf file to search.
#' @param diag_ion The mz of the diagnsotic ion(s)
#' @param tol A mass tolerance to respect during the search
#' @param export_mgf \code{logical} Export only scans containing diagnsotic ion(s) into a an mgf file.
#' @param save_file \code{logical} Save the results as csv file
#' @return A \code{tibble} with 6 columns including the diagnostic ion \code{m/z} and its intensity relative to the base peak.
#' @importFrom dplyr near
#' @importFrom tibble tibble
#' @import Spectra
#' @importFrom MsBackendMgf MsBackendMgf
#' @importFrom cli cli_abort cli_progress_bar cli_progress_update cli_progress_done cli_inform cli_alert_info cli_alert_success
#' @importFrom rlang is_missing
#' @importFrom BiocManager install
#' @importFrom utils install.packages write.csv
#'
#' @export
mgf_searchDiagIon <- function(mgf_file, diag_ion, tol = 0.002, export_mgf = FALSE, save_file = FALSE){



# check inputs --------------------------------------------
  if(rlang::is_missing(mgf_file)){
    cli::cli_abort("No mgf file name is identified.")
  }
  if(rlang::is_missing(diag_ion)){
    cli::cli_abort("No diagnostic ion is identified.")
  }

  #to avoid any non-numeric input
  stopifnot( "'diag_ion' must be a number or a numeric vector." = is.numeric(diag_ion) == TRUE)

####

# check dependencies availability -----------------------

if (!requireNamespace("BiocManager")) install.packages("BiocManager")

  if (!requireNamespace("Spectra", quietly = TRUE))
    BiocManager::install("Spectra")

  if (!requireNamespace("MsBackendMgf", quietly = TRUE))
    BiocManager::install("MsBackendMgf")

  if (!requireNamespace("BiocParallel", quietly = TRUE))
    BiocManager::install("BiocParallel")



####---


# iterate over mgf files --------------------
for (file in seq_along(mgf_file)){

file_name = tools::file_path_sans_ext(basename(mgf_file))

## read mgf -----------------
  sps <- .mgf_to_sp(mgf_file[file])


  msg = paste0("Extracting diagnostic ions from ", basename(mgf_file[file]))
  cli::cli_progress_bar(type = "iterator", name = msg )

  num_spectra <- length(Spectra::acquisitionNum(sps))

  # Pre-allocate matrices: one for double values and one for integer values
  double_matrix <- matrix(NA_real_, nrow = 0, ncol = 4)  # diag_ion, diag_relint, prec_mz, rt
  int_matrix <- matrix(NA_integer_, nrow = 0, ncol = 2)  # scan, prec_z

  colnames(double_matrix) <- c("diag_ion", "diag_relint", "prec_mz", "rt")
  colnames(int_matrix) <- c("scan", "prec_z")

  # Pre-allocate vector for file names (character)
  file_vec <- character()

### iterate over spectra in mgf -----------------

  for (i in seq_len(num_spectra)) {

    temp_double <- matrix(NA_real_, nrow = 0, ncol = 4)  # temp matrix for double values
    temp_int <- matrix(NA_integer_, nrow = 0, ncol = 2)  # temp matrix for integer values

#### iterate over user-defined diag ions -----------------
     for (ion in seq_along(diag_ion)) {

      isIonhere <- dplyr::near(diag_ion[ion], Spectra::mz(sps)[[i]], tol = tol) #logical vector

##### iterate over same spectrum containing diag ions falling within tolerance range -----------------

          #find at least one is TRUE
          if (any(isIonhere, na.rm = TRUE)) {
            indices =   which(isIonhere) #find its/their index/indices

            # Allocate temp matrices for this spectrum/ion
            temp_double_ion <- matrix(NA_real_, nrow = length(indices), ncol = 4)
            temp_int_ion <- matrix(NA_integer_, nrow = length(indices), ncol = 2)

              for (indx in seq_along(indices)) {
###### fill in data with some manipulations -----------------

                  temp_double_ion[indx, ] <- c(
                    round(Spectra::mz(sps)[[i]][indices[indx]], 4),
                    round(Spectra::intensity(sps)[[i]][indices[indx]] / max(Spectra::intensity(sps)[[i]]),3),
                    round(Spectra::precursorMz(sps)[[i]], 4),
                    round(Spectra::rtime(sps)[[i]] / 60, 2) #sec to min
                  )

                  temp_int_ion[indx, ] <- c(
                    .found_replace_scan(sps, i),  # scan number
                    Spectra::precursorCharge(sps)[[i]]  # precursor charge
                  )

              }

            # Bind rows to temp matrices
            temp_double <- rbind(temp_double, temp_double_ion)
            temp_int <- rbind(temp_int, temp_int_ion)
          }
     }

    # Append the current spectrum's diagnostic ions to the main matrices
    double_matrix <- rbind(double_matrix, temp_double)
    int_matrix <- rbind(int_matrix, temp_int)

    # Add the file name corresponding to this spectrum
    if (nrow(temp_double) > 0) {
      file_vec <- c(file_vec, rep(file_name, nrow(temp_double)))
    }

    cli::cli_progress_update()
  }

  if(nrow(double_matrix) == 0L){cli::cli_alert_info("No diagnostic ion was found.")}else{

    # Combine the matrices and vector into a data frame
    final_df <- tibble::tibble(
      file = file_vec,
      diag_ion = double_matrix[, "diag_ion"],
      diag_relint = double_matrix[, "diag_relint"],
      prec_mz = double_matrix[, "prec_mz"],
      rt = double_matrix[, "rt"],
      scan = int_matrix[, "scan"],
      prec_z = int_matrix[, "prec_z"],
      stringsAsFactors = FALSE
    )

    cli::cli_progress_done()
}



#no need to save a file if no diagnostic ion was found.
  if(save_file && nrow(final_df) != 0L){
    diagnostic_ions <- paste(as.integer(diag_ion), collapse = "_")
    file_csv = paste0("diagIons_",diagnostic_ions, "_", file_name, ".csv")
    write.csv(x = final_df, file = file_csv, row.names = FALSE)
    wd = getwd()
    cli::cli_alert_success("'{file_csv}' is saved sucessfully into '{wd}' directory.")
  }


  # Export scans as containing specific diag ions ----------------------

  if (export_mgf)
  #no scans at all

     if(nrow(final_df) == 0L) {
    cli::cli_alert_info("No scan number was found.")
  } else {

    # scans are all NA
       if (all(is.na(sps$acquisitionNum))) {
           for (i in seq_along(sps)) {
        sps$acquisitionNum[i] = .found_replace_scan(sps, i)
      }
    }

    diagnostic_ions <- paste(as.integer(diag_ion), collapse = "_")
    file_submgf <- paste0("selectedscans_", diagnostic_ions, "_", file_name, ".mgf")
    sub_scans <- final_df$scan
    sub_sps <- sps[sps$acquisitionNum %in% sub_scans]
    map <- Spectra::spectraVariableMapping(MsBackendMgf::MsBackendMgf())

    Spectra::export(
      sub_sps,
      backend = MsBackendMgf::MsBackendMgf(),
      file = file_submgf,
      mapping = map
    )

    cli::cli_alert_success("'{file_submgf}' is saved successfully.")
  }




  }






#add a space between messages and output df
cat("\n")

if(nrow( final_df) != 0L) {return(final_df)

}
invisible(NULL)  # Return NULL if no results
}



#' @examples
#' # example code
#'



#' @noRd
.mgf_to_sp <- function(mgf_file){

  mgf_file_basename = basename(mgf_file)

  cli::cli_inform(message = paste0("Reading ", mgf_file_basename, "..."))
  BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
  mgftospec <- suppressMessages(Spectra::Spectra(object = {{mgf_file}},
                                                 source = MsBackendMgf::MsBackendMgf()))

  cli::cli_alert_success( "'{mgf_file_basename}' is converted into a 'Spectra' object.")

return(mgftospec)
}



#' @noRd

.found_replace_scan = function(spec_obj, scan_i){

  ifelse(!is.na(spec_obj$acquisitionNum[scan_i]),
         spec_obj$acquisitionNum[scan_i], #return the scan
         {
           #otherwise extract the scan number from the title
           title = spec_obj[['TITLE']][scan_i]

           pattern <- "(scan=|_?scan:|Scan|Index:)\\s*(\\d+)"

           # Search for the pattern in the title
           match <- as.integer(stringr::str_extract(title, pattern = pattern, group = 2))


           spec_obj$acquisitionNum[scan_i] <- match

          #if match is integer that's fine, any other type will be coerced to NA.
           return(match)
         }
  )
}

