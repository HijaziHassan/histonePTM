




.normalize_vec <- function(x){

  norm_x <-  x /sum(x, na.rm = TRUE)

  return(norm_x)
}



#' A function to normalize intensities of peptidoforms by dividing an intensity value by the
#' sum of all intensities other peptidoforms of the same sequence.
#'
#' @param df A \code{tibble} or a \code{dataframe}
#' @param select_cols Either one of the select helper functions with proper argument
#' (e.g. \code{starts_with("WT_")}) or a vector \code{c()}containing the names of
#'  intensity/abundance columns to be normalized.
#'
#' @return The provided \code{dataframe} is returned with relative instead of absolute
#' intensities values of the abundance/intensity columns passed to the function.
#'
#' \deqn{\text{Relative Intensity} = \frac{\text{Intensity}_{\text{peptide X}}}{\sum \text{Intensity}_{\text{all peptides of the same sequence as X}}}}
#'
#'
#' @examples
#' quant_relIntensity(mtcars, select_cols = contains("m"))
#' quant_relIntensity(mtcars, select_cols = c(wt,  qsec, vs))
#'
#'@import rlang
#'@import dplyr
#'
#' @export
quant_relIntensity <- function(df, select_cols){

  select_exp <- rlang::enexpr(select_cols)

  df_norm <- df |>
    dplyr::mutate( dplyr::across(!!select_exp, ~ .normalize_vec(x=.))
    )



  rlang::eval_tidy(df_norm)

}











