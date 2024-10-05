
#' Rename 'abundance_xxx' columns in Proline Export.
#' @description
#'  Create custom names to replace dat file coded 'abundance' and 'psm count' headers in Proline export.
#'
#' @param raw_df Dataframe extracted using \code{misc_extractMetaData} function.
#' @param user_df Dataframe containing "SampleName" and file" columns (could also contain 'Condition', 'Bioreplicate' or 'TechReplicate' columns).
#' @importFrom dplyr left_join mutate any_of select bind_rows
#' @importFrom tibble deframe
#' @return named vector
#' @export
#'
misc_useCustomNames <- function(raw_df, user_df){


  #remove unnammed and automatically renamed columns "ex: ...1" before merging.
  raw_df <- raw_df[, !grepl("^\\.{3}[1-9]$|^$", colnames(raw_df))]
  user_df <- user_df[, !grepl("^\\.{3}[1-9]$|^$", colnames(user_df))]


  # the unmatced argument is not working.
  abn_ColNames <- dplyr::left_join(raw_df,
                                   user_df,
                                   by= "file",
                                   keep = FALSE,
                                   unmatched = "drop") |>
    dplyr::mutate(new_name = paste0("abundance_", channel),
           .keep = "unused")


  psmcount_ColNames <- dplyr::left_join(raw_df,
                                 user_df ,
                                 by= "file",
                                 keep = FALSE,
                                 unmatched = "drop") |>
    dplyr::mutate(SampleName = paste0(SampleName, "_psmCount"),
           new_name = paste0("psm_count_", channel),
           .keep = "unused")

  ColNames <- dplyr::bind_rows(abn_ColNames, psmcount_ColNames) |>
    dplyr::select(-c(file, dat, dplyr::any_of(c('Condition', 'condition', 'Bioreplicate', 'BioReplicate', 'TechReplicate')))) %>%
    tibble::deframe()

  return(ColNames)


}
