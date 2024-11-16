
#' Extract raw file and search file names
#' @description
#' Extract raw and search result files names from Proline output excel file (\code{"Search settings and infos"} sheet).
#'
#' @param filename Proline export excel filename ('\code{Search settings and infos}' sheet is hardcoded.)
#' @param save_file Provide a file name if you want to save the extracted metadata.
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter rename mutate
#' @importFrom stringr str_remove
#' @importFrom readxl read_excel
#' @importFrom openxlsx write.xlsx
#' @importFrom utils write.csv
#'
#' @return \code{.csv} file with 3 columns: \code{file}, \code{channel}, and \code{dat}.
#' @export
misc_extractMetaData <- function(filename, save_file= ""){

  suppressWarnings({

    extract_names <- readxl::read_excel(path = {{filename}},
                                        sheet = "Search settings and infos",
                                        .name_repair = "unique_quiet") |>
      #I used those row identifiers since they'are always present even if one sample is being processed in Proline.
      dplyr::filter(project_name %in% c("result_file_name", "raw_file_name", "job_number")) |>
      tibble::as_tibble()

  })

  dat_files <- tibble::as_tibble(t(extract_names[ , - 1]), .name_repair = "minimal")

  colnames(dat_files) = c('file', 'channel', 'dat')

  dat_files <- dat_files |>
    dplyr::mutate(dat= as.integer(dat)) |>
    dplyr::mutate(channel= stringr::str_remove(channel, ".dat|.pep.xml|.pepXML|.msf|.pdResult|.mzid|.msr"))

  if(save_file != ""){
    if(tools::file_ext(save_file) %in% c('xls', 'xlsx')){

    openxlsx::write.xlsx({{save_file}})
    }else{
      if(tools::file_ext(save_file) == ""){

        save_file = paste0(save_file, ".csv")
      }
  write.csv(x = dat_files, file = {{save_file}}, row.names = FALSE)
      }
    }

  return(dat_files)
}


