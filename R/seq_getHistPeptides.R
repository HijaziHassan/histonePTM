



#' Extract histone peptides
#' @description
#' Isolate specific fragments of user-defined histone protein(s).
#'
#' @param df A \code{dataframe} containing the proteomics experiment results.
#' @param seq_col The name of the column containing peptide sequences.
#' @param histoneProtein Histone protein(s) to be studied: "H3", "H4" or both "H3-H4".
#'
#' @return A subset of \code{df} only containing peptides of the defined histone protein(s) with iRTs if present.
#'
#' @import rlang
#' @import purrr
#' @import dplyr
#' @export

seq_getHistPeptide <- function(df, seq_col = sequence ,  histoneProtein = c("H3-H4", "H3", "H4")){

  seq_col = rlang::enexpr(seq_col)
  #seq_col argument must be unquoted
  seq_col <-  ifelse(is.character(seq_col),
          rlang::ensym(seq_col),
         rlang::enexpr(seq_col)
         )



  histoneProtein = match.arg(histoneProtein)

  H3 <- purrr::pluck(histonePTM::sequenceDB, "H3", "argC_frags")
  H4 <- purrr::pluck(histonePTM::sequenceDB, "H4", "argC_frags")
  iRT <- purrr::pluck(histonePTM::sequenceDB, "iRT", "tryp_frags")

  if(histoneProtein == "H3"){
    df <- df |>
      dplyr::filter({{seq_col}} %in% c(H3,iRT)) #iRT added
  }

  else if (histoneProtein == "H4") {

    df <- df |>
      dplyr::filter({{seq_col}} %in% c(H4, iRT)) #iRT added

  }else{

    df <- df |>
      dplyr::filter({{seq_col}} %in% c(H3, H4, iRT))


  }


  return(df)
}





