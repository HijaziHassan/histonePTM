#' Extract user-defined Scan from MGF file.
#' @description
#' wrapper function around \code{mgf.read} from \code{spectrum_utils} library in Python.
#'
#' @param mgf_file mgf file (path) name
#' @param scan \code{numeric}. scan number
#'
#' @return tibble with 4 colums: file name, scan number and their corresponfing mz and intensity values.
#' @export

mgf_extractScan <- function(mgf_file, scan){

reticulate::source_python("./inst/python/extract_scan_mgf.py")

  spec <- extract_scan(mgf_file_path = mgf_file, target_scan = scan)

  if(!is.null(spec)){


  file_name = basename(mgf_file)

  my_spec = data.frame(
    mgf = file_name,
    scan = spec$params$scans,
    mz = spec$`m/z array`,
    intenisty = spec$`intensity array`
  ) |> tibble::as_tibble()

return(my_spec)
}
}

