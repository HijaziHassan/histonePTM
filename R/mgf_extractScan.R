#' Extract user-defined Scan from MGF file.
#' @description
#' wrapper function around \code{mgf.read} from \code{spectrum_utils} library in Python.
#'
#' @param mgf_file mgf file (path) name
#' @param scan \code{numeric}. scan number
#'
#' @importFrom dplyr bind_rows
#' @importFrom tibble as_tibble
#' @importFrom cli cli_abort
#'
#' @return tibble with 4 colums: file name, scan number and their corresponfing mz and intensity values.
#' @export

mgf_extractScan <- function(mgf_file, scan){

  if (!file.exists(mgf_file)) {
    cli::cli_abort("The file {.arg {mgf_file}} was not found.")
  }


  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  My_Package_Name_path <- system.file("python", package = "histonePTM")
  reticulate::py_run_string(paste0("import sys; sys.path.append('", My_Package_Name_path, "')"))

  reticulate::source_python(system.file("python",
                                        "extract_scan_mgf.py",
                                        package = "histonePTM",
                                        mustWork = TRUE
  ))


#reticulate::source_python("./inst/python/extract_scan_mgf.py")

list_df <- vector(mode = "list", length = length(mgf_file)*length(scan))


for(mgf in mgf_file){

  file_name = basename(mgf_file)

  for(s in scan){

spec <- extract_scan(mgf_file_path = mgf, target_scan = s)

  if(!is.null(spec)){

    my_spec = data.frame(
    mgf = file_name,
    scan = as.integer(spec$params$scans),
    mz = spec$`m/z array`,
    intensity = spec$`intensity array`
              ) |>
      tibble::as_tibble()


  # Append the data frame to the results list
  list_df[[length(list_df) + 1]] <- my_spec

          }
      }

}

  final_df <- dplyr::bind_rows(list_df)

  return(final_df)


  }
