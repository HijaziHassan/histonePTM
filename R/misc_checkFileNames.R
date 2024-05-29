# Check file names provided with those actually found in the Proline excel output.
misc_checkFileNames <- function(raw_df, user_df, common_col) {
  check.presence = all(raw_df[[common_col]] %in% user_df[[common_col]])
  check.extras = all.equal(sort(raw_df[[common_col]]), sort(user_df[[common_col]]))
  diff_col = setdiff(user_df[[common_col]], raw_df[[common_col]])

  if(isFALSE(check.presence)){
    cli::cli_abort(c("Mismatch",
                     "i"= "The provided name(s) must match the one(s) presented in 'Search settings and infos' sheet." ,
                     "x"= "You either did not add or misspelled  the following filename(s): {diff_col}."))
  }

  if(isTRUE(check.extras)){
    cli::cli_alert_success(cat("\nThe provided and actual names match perfectly. Good job!"))
  }else{cli::cli_alert_warning(cat("\nExtra unmatched names are dropped"))}
}
