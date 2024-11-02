
#' @title Normalize Intensity values
#'
#' @description A function to normalize intensities of peptidoforms by dividing an intensity value by the
#' sum of all intensities other peptidoforms of the same sequence.
#'
#' @param df A \code{tibble} or a \code{dataframe}
#' @param select_cols Either one of the select helper functions with proper argument
#' (e.g. \code{starts_with("WT_")}) or a vector \code{c()}containing the names of
#'  intensity/abundance columns to be normalized.
#' @param grouping_var A unqiue variable to group by like 'peptide sequence' or so. This mandatory if `norm_method` is 'peptide_family".
#' @param norm_method Normalization method. Either by peptide family or by total intensity of all (histone) peptides in the sample. The latter depends on what is in your dataset and if you have prefiltered it or not.
#'
#' @return The provided \code{dataframe} is returned with relative instead of absolute
#' intensity values of the abundance/intensity columns passed to the function.
#'
#' \deqn{\text{norm_peptide_family} = \frac{\text{Intensity}_{\text{peptide X}}}{\sum \text{Intensity}_{\text{all peptides of the same sequence as X}}}}
#' \deqn{\text{norm_peptide_family} = \frac{\text{Intensity}_{\text{peptide X}}}{\sum \text{Intensity}_{\text{all peptides in the sample}}}}
#'
#'
#' @examples
#' quant_relIntensity(iris, contains("Length"), grouping_var = "Species")
#' quant_relIntensity(iris, starts_with("Sepal"), grouping_var = "Species")
#'
#' @importFrom dplyr across mutate
#'
#' @export

quant_relIntensity <- function(df, select_cols, grouping_var, norm_method = c('peptide_family', "total_histone")){

  norm_by = match.arg(norm_method)


  if ( norm_by == 'peptide_family'){

  df_norm <- df |>
    dplyr::mutate(dplyr::across({{select_cols}}, ~ .normalize_vec(x=.)),
                   .by = {{grouping_var}}
    )

  }else if(norm_by =="total_histone" ){

    df_norm <- df |>
      dplyr::mutate(dplyr::across({{select_cols}}, ~ .normalize_vec(x=.)))

}

  return(df_norm)

}


#' @noRd

.normalize_vec <- function(x){

  norm_x <-  x /sum(x, na.rm = TRUE)

  return(norm_x)
}






