
#' Extract raw file and search file names
#' @description
#' Extract .raw file and .dat file names from Proline output excel file (\code{"Search settings and infos"} sheet).
#'
#' @param filename Proteomics experiment excel output.
#' @import tibble
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import readxl
#'
#' @return \code{.csv} file with 3 columns: \code{file}, \code{channel}, and \code{dat}.
#' @export
extract_metadata <- function(filename, save_file= FALSE){

  suppressWarnings({

    extract_names <- readxl::read_excel(path = {{filename}},
                                        sheet = "Search settings and infos",
                                        .name_repair = "unique_quiet") %>%
      #I used those row identifiers since they'are always present even if one sample is being processed in Proline.
      dplyr::filter(project_name %in% c("result_file_name", "raw_file_name", "job_number")) %>%
      tibble::as_tibble()

  })

  dat_files <- tibble::as_tibble(t(extract_names[ , - 1])) %>%
    dplyr::rename(file = V1, channel = V2, dat = V3) %>%
    dplyr::mutate(dat= as.integer(dat)) |>
    dplyr::mutate(channel= stringr::str_remove(channel, ".dat"))

  if(save_file != ""){
  write.csv(x = dat_files, file = {{save_file}}, row.names = FALSE)
  }

  return(dat_files)
}


