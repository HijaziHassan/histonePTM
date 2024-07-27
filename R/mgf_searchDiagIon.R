
#' Find MS/MS spectra that contain a diagnostic ion.
#' @description
#' Report user-defined diagnsotic ion(s) per MS/MS spectrum from an \code{mgf} file.
#' @param mgf_file mgf file to search.
#' @param diag_ion The mz of the diagnsotic ion(s)
#' @param tol A mass tolerance to respect during the search
#' @param export_mgf \code{logical} Export only scans containing diagnsotic ion(s) into a an mgf file.
#' @param save_file \code{logical} Save the results as csv file
#'
#' @return A \code{tibble} with 6 columns including the diagnostic ion \code{m/z} and its intensity relative to the base peak.
#' @import dplyr
#' @import Spectra
#' @importFrom MsBackendMgf MsBackendMgf
#' @import Spectra
#' @import cli
#' @import rlang
#' @importFrom BiocManager install
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
  sapply(libraries, requireNamespace, character=TRUE)

####---


# iterate over mgf files --------------------
for (file in seq_along(mgf_file)){

file_name = tools::file_path_sans_ext(basename(mgf_file))

## read mgf -----------------
  sps <- .mgf_to_sp(mgf_file)


  msg = paste0("Extracting diagnostic ions from ", basename(mgf_file))
  cli::cli_progress_bar(type = "iterator", name = msg )

  final_list = list()

### iterate over spectra in mgf -----------------

  for (i in 1:length(Spectra::acquisitionNum(sps))) {

    temp_all_diag <- vector(mode = 'list')

#### iterate over user-defined diag ions -----------------
     for (ion in seq_along(diag_ion)) {

      isIonhere <- dplyr::near(diag_ion[ion], Spectra::mz(sps)[[i]], tol = tol) #logical vector

##### iterate over same spectrum containing diag ions falling within tolerance range -----------------

          #find at least one is TRUE
          if (any(isIonhere, na.rm = TRUE)) {
            indices =   which(isIonhere) #find its/their index/indices

            temp_one_diag <- list() # to account for two values too close to each other.

              for (indx in seq_along(indices)) {
###### fill in data with some manipulations -----------------

                  temp_one_diag[[indx]] <- list(
                    file         = file_name,
                    diag_ion     = round(Spectra::mz(sps)[[i]][indices[indx]], 4),
                    diag_relint  = round(Spectra::intensity(sps)[[i]][indices[indx]] / max(Spectra::intensity(sps)[[i]]),3),
                    #scan         = Spectra::acquisitionNum(sps)[[i]],
                    scan         = .found_replace_scan(sps, i),
                    prec_mz      = round(Spectra::precursorMz(sps)[[i]], 4),
                    prec_z       = Spectra::precursorCharge(sps)[[i]],
                    rt           = round(Spectra::rtime(sps)[[i]] / 60, 2) #sec to min
                  )




              }

            cli::cli_progress_update()

        # Update the combined list
        temp_all_diag[[ion]] <- temp_one_diag


      }


    }


    final_list[[i]] <- dplyr::bind_rows(temp_all_diag)


  }


  final_df <- dplyr::bind_rows(final_list)

  if(nrow( final_df) == 0L){cli::cli_alert_info("No diagnostic ion was found.")}else{

    final_df <- final_df |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("diag_ion", "diag_relint", "prec_mz", "rt")), as.double),
      dplyr::across(dplyr::all_of(c("scan", "prec_z")), as.integer)
    )

    cli::cli_progress_done()
}



#no need to save a file if no diagnsitic ion was found.
  if(save_file && nrow(final_df) != 0L){
    diagnostic_ions <- paste(as.integer(diag_ion), collapse = "_")
    file_csv = paste0("diagIons_",diagnostic_ions, "_", file_name, ".csv")
    write.csv(x = final_df, file = file_csv, row.names = FALSE)
    wd = getwd()
    cli::cli_alert_success("'{file_csv}' is saved sucessfully into '{wd}' directory.")
  }


  # Export scans as containing specific scans ----------------------

    if (export_mgf)
  #no scans at all

    if (nrow(final_df) == 0L) {
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

if(nrow( final_df) != 0L) return(final_df)

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

  dplyr::if_else(!is.na(spec_obj$acquisitionNum[scan_i]),
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

