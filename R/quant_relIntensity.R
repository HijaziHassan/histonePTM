




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
#' quant_relIntensity(iris, contains("Length"), grouping_var = "Species")
#' quant_relIntensity(iris, starts_with("Sepal"), grouping_var = "Species")
#'
#'@import rlang
#'@import dplyr
#'
#' @export
# quant_relIntensity <- function(df, select_cols, grouping_var){
#
#   select_exp <- rlang::enexpr(select_cols)
#
#   df_norm <- df |>
#     dplyr::mutate( dplyr::across(!!select_exp, ~ .normalize_vec(x=.)),
#                    .by = {{grouping_var}}
#     )
#
#
#
#   rlang::eval_tidy(df_norm)
#
# }
quant_relIntensity <- function(df, select_cols, grouping_var){

  #select_exp <- rlang::enexpr(select_cols)

  df_norm <- df |>
    dplyr::mutate( dplyr::across({{select_cols}}, ~ .normalize_vec(x=.)),
                   .by = {{grouping_var}}
    )



  return(df_norm)

}


quant_relIntensity(iris, contains("Length"), grouping_var = "Species")








