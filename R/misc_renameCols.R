
#' Rename columns in Proline Output
#' @description
#' Rename columns "abundance_..." from Proline proteomics results output.
#'
#'
#' @param analysisfile Name of Proline output excel file
#' @param metafile Name of user-defined sample names (at least two columns: SampleName [user defined sample names] and file [raw files]).
#'
#' @return renamed dataframe
#' @importFrom dplyr rename any_of
#' @importFrom readxl read_xlsx
#' @export
#'
#'
misc_renameCols <- function(analysisfile, metafile){


  df_raw <- readxl::read_xlsx(analysisfile, sheet= 'Best PSM from protein sets')
  meta_raw <-  misc_extractMetaData(analysisfile)
  meta_user <- readxl::read_xlsx(metafile)

  #combine
  names_vector <- misc_useCustomNames(raw_df = meta_raw, user_df = meta_user)

  #compare
  misc_checkFileNames(meta_raw, user_df = meta_user, common_col = "file")

  #rename

  renamed_df <- df_raw |> dplyr::rename(!!!dplyr::any_of(names_vector))


  return(renamed_df)

}


