



#' Extract selected histone peptides
#' @description
#' Isolate specific fragments of user-defined histone protein(s).
#'
#' @param df A \code{dataframe} containing the proteomics experiment results.
#' @param seq_col The name of the column containing peptide sequences.
#' @param histoneProtein Histone protein(s) to be studied: "H3", "H4" or both "H3-H4".
#'
#' @return A subset of \code{df} only containing peptides of the defined histone protein(s) with \href{https://biognosys.com/product/irt-kit/}{Biognosys}
#' iRTs if present.
#'
#' @importFrom rlang ensym enexpr
#' @importFrom purrr pluck map
#' @importFrom dplyr filter
#' @importFrom cli cli_abort
#' @export

seq_getHistPeptide <- function(df, seq_col = sequence ,  histoneProtein = c("All", "H3", "H4", "H2A", "H2B", "H1")){

  seq_col = rlang::enexpr(seq_col)
  #seq_col argument must be unquoted
  seq_col <-  ifelse(is.character(seq_col),
          rlang::ensym(seq_col),
         rlang::enexpr(seq_col)
         )


  if(!any(histoneProtein %in% c("All", "H3", "H4", "H2A", "H2B", "H1"))){cli::cli_abort(c("Invalid input.",
                                                                               "x" = "'{histoneProtein}' is not recognized.",
                                                                               "i" = "Please choose any of: 'All' (default), 'H3', 'H4', 'H2A', 'H2B', 'H1')"))
 }else if ("All" %in% histoneProtein){

   histoneProtein = c("H3", "H4", "H2A", "H2B", "H1")
   protein_extracted <- purrr::map(.x = histoneProtein,
                                   .f = ~ purrr::pluck(histonePTM::sequenceDB, .x, "argC_frags"))|>
                                         unlist()

    iRT <- purrr::pluck(histonePTM::sequenceDB, "iRT", "tryp_frags")

    df <- df |>
      dplyr::filter({{seq_col}} %in% c(protein_extracted, iRT))


 }else{


   protein_extracted <- purrr::map(.x = histoneProtein,
                                   .f = ~ purrr::pluck(histonePTM::sequenceDB, .x, "argC_frags")) |>
                                       unlist()

   iRT <- purrr::pluck(histonePTM::sequenceDB, "iRT", "tryp_frags")

   df <- df |>
     dplyr::filter({{seq_col}} %in% c(protein_extracted, iRT))

 }

  cli::cli_inform("{histoneProtein} histone peptides were extracted.")

  return(df)
}






