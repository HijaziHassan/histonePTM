
#' Site-specific abundance
#' Calculates the abundance of an individual histone PTM site.
#' @param df Dataframe containing normalized values
#' @param df_meta Dataframe containing the sample design. Must contain the column "SampleName" and "Condition".
#' @param ptm_col PTM column containing PTMs in Brno nomenclature, separated by "-" if multiple PTMs exist (e.g. 'K27ac-K36bu')
#' @param id_col The identity column(s) (e.g. column containing sequences).
#' @param int_cols Column(s) containing intensity values.
#' @param format Format of `df` ('wide' (default) or 'long'))
#' @param remove_ptm A ptm to be excluded. provided as a string.
#' @param save_plot \code{bool, FALSE} (default)
#'
#' @importFrom dplyr select filter pull summarize across bind_rows
#' @importFrom stringr str_detect
#' @importFrom purrr map
#' @return dataframe with id_col(s), `PTMsite` containing the quantified sites, and the intensity column(s)
#' @export
#'

ptm_siteAbundance<- function(df, metadata, ptm_col, id_col, int_cols, format= c('wide', 'long'), remove_ptm = NULL, save_plot = FALSE) {

  format= match.arg(format)

  #for removing Nterm 'prNt-' or '*' from flagging duplications function
  if (is.null(remove_ptm)) {
    remove_ptm <- "\\*"
  } else {
    remove_ptm <- paste("\\*", remove_ptm, sep = "|")
  }

  #extact unique residue and ptm (e.g. K10ac)
  res_ptm <- df |>
    dplyr::pull({{ ptm_col }}) |>
    strsplit(split = "-") |>
    unlist() |>
    unique() |>
    gsub(pattern = remove_ptm, replacement = "", x = _) |>
    (\(x) x[x != "" & !is.na(x)])()



  #iterate over unique ptm list, sum intensities
  df_siteAbundance <- purrr::map(.x = res_ptm, .f = ~ {
    df |>
      dplyr::filter(stringr::str_detect({{ptm_col}}, .x)) |>
      dplyr::summarize(
        PTMsite = .x,
        dplyr::across({{int_cols}}, .fns= \(ptm) sum(ptm, na.rm = TRUE)),
        .by = {{id_col}}
      )
  }) |> dplyr::bind_rows()





  if(save_plot){

    if(format == "wide"){

      df_plot <- df_siteAbundance |>
      tidyr::pivot_longer(cols = -c(rlang::ensym(id_col), PTMsite),
                          names_to= 'SampleName',
                          values_to = 'intensity') |>
      dplyr::left_join(y = df_meta, by = 'SampleName')
}

   plot_jitterbarIntvsPTM(dataset = df_plot
                                   , x_axis = PTMsite
                                   , y_axis = intensity
                                   , condition = Condition
                                   , id_col = {{id_col}}
                                   , save_plot = save_plot)
}


  return(df_siteAbundance)
}

