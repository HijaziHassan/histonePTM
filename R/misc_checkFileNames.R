#' Check match betwen two columns
#' @description
#' Could be used to check for user input file names if they match file names extracted from 'Search settings and infos' sheet of Proline export.
#'
#' @param raw_df Dataframe containing acolumn of correct file names
#' @param user_df Dataframe containing a column of file names provided by the user to which another custom name are assigned
#' @param common_col column of file names in both dataframes (should be the same name in both).
#'
#' @return error or success message.
#' @export
#' @examples
#' df1 <- data.frame(file= c('HF1_022210', 'HF1_022222'))
#' df2 <- df1
#' misc_checkFileNames(df1, df2, common_col = 'file')
#' #>The provided and actual names match perfectly. Good job!✔



misc_checkFileNames <- function(raw_df, user_df, common_col) {

  if(missing(common_col)) stop('"common col" argument must be provided.')

  #remove columns with no names or column renamed after repair to ...1, ..2
  raw_df <- raw_df[, !grepl("^\\.{3}[1-9]$|^$", colnames(raw_df))]
  user_df <- user_df[, !grepl("^\\.{3}[1-9]$|^$", colnames(user_df))]

  check.presence = all(raw_df[[common_col]] %in% user_df[[common_col]])
  check.extras = all.equal(sort(raw_df[[common_col]]), sort(user_df[[common_col]]))
  diff_col = setdiff(user_df[[common_col]], raw_df[[common_col]])

  if(isFALSE(check.presence)){
    cli::cli_abort(c("Mismatch",
                     "i"= "The provided name(s) must match the one(s) presented in 'Search settings and infos' sheet." ,
                     "x"= "You either did not add or misspelled  the following filename(s): {diff_col}."))
  }

  if(isTRUE(check.extras)){
    cli::cli_alert_success("The provided and actual names match perfectly. Good job!")
  }else{cli::cli_alert_warning("There are extra unmatched names")

  return(diff_col)
    }

}
