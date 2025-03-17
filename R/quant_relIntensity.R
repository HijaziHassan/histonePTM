
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
#' @param norm_method Normalization method. Either by 'peptide_family' or by 'peptide_total' intensity. The latter depends on what is in your dataset and if you have prefiltered it or not.
#'
#' @return The provided \code{tibble} or \code{dataframe} is returned with relative instead of absolute
#' intensity values of the abundance/intensity columns passed to the function.
#'
#' \deqn{ norm\_peptide\_family = \frac{Intensity_{peptide\ X}}{\sum Intensity_{all\ peptides\ of\ the\ same\ sequence\ as\ X}} }{
#' norm_peptide_family = Intensity_{peptide X} / sum(Intensity_{all peptides of the same sequence as X}) }
#'
#' \deqn{ norm\_peptide\_total = \frac{Intensity_{peptide\ X}}{\sum Intensity_{all\ peptides\ of\ the\ same\ sequence\ as\ X}} }{
#' norm_peptide_family = Intensity_{peptide X} / sum(Intensity_{all\ peptides\ of\ all\ filtered\ peptdies\ in\ each\ sample}) }

#'
#'
#' @examples
#' \dontrun{
#' #dummy examples
#' quant_relIntensity(iris, contains("Length"), grouping_var = "Species")
#' quant_relIntensity(iris, starts_with("Sepal"), grouping_var = "Species")
#'}
#' @importFrom dplyr across mutate all_of
#'
#' @export

quant_relIntensity <- function(df, select_cols, grouping_var, norm_method = c('peptide_family', "peptide_total")){

  norm_by = match.arg(norm_method)


  if ( norm_by == 'peptide_family'){

  df_norm <- df |>
    dplyr::mutate(dplyr::across(dplyr::all_of({{select_cols}}), ~ .normalize_vec(x=.)),
                   .by = {{grouping_var}}
    )

  }else if(norm_by =="peptide_total" ){

    df_norm <- df |>
      dplyr::mutate(dplyr::across(dplyr::all_of({{select_cols}}), ~ .normalize_vec(x=.)))

}

  return(df_norm)

}


#' @noRd

.normalize_vec <- function(x){

  norm_x <-  x /sum(x, na.rm = TRUE)

  return(norm_x)
}






