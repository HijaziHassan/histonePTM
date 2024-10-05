
#' Assessing Monomethyl Derivatization
#'
#' @description
#' Evaluate the efficiency of monomethylated lysine derivarization per sequencexPTM combination.
#'
#' @param df Dataframe containing proteomic analysis dataset
#' @param seq_col Sequence column
#' @param ptm_col PTM column
#' @param int_col Intensity column(s)
#' @param cond_col (`optional`) Condition column (e.g. concentration) used for legend. If absent
#' @param me1 Monomethyl representation inside the `ptm_col`column (i.e. 'me1' (`default`) or 'methyl (K)', ...)
#' @param me1_label Monomethyl+label representation inside the `ptm_col`column (i.e. 'bu' (`default`) or 'Butyryl (K)', ....)
#' @param format if set to `wide` (most often), data will be reshaped into long format before plotting. Check `df_meta`.
#' @param df_meta (`optional`). If df is in `wide` format, then one can provide df_meta containing at least two columns (`SampleNames` and `Conditions`).
#' `SampleNames` MUST match exactly the names of the `int_col`columns.
#' @param save_plot (`logical`). If `TRUE`, all the plots will be saved in one pdf. If `FALSE` (`default`), the plot(s) will be shown on Plots pane.)
#'
#' @importFrom dplyr select left_join filter mutate inner_join right_join pull bind_rows
#' @importFrom stringr str_detect str_replace_all
#' @importFrom cli cli_alert_warning cli_alert_warning col_red cli_abort
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom grDevices dev.off pdf
#' @return Boxplot with jittered point to assess me1 labeling when using anhydrides to label lysines.
#' @export
#'

quant_me1LabelingPlot <- function(df,
                                  seq_col,
                                  ptm_col,
                                  int_col,
                                  cond_col,
                                  me1 = "me1",
                                  me1_label = "bu",
                                  format="",
                                  df_meta,
                                  save_plot = FALSE){



  #Avoid forgetting or using using wrong pattern inside tidyselect function

 df_check <- df |> dplyr::select({{int_col}})

 if (ncol(df_check) == 0) {cli::cli_abort(c("Columns not found!",
   'x' = "The provided column names in `int_col` argument are not found in `df`.",
                                            "i" = "provide correct names using `tidyselect` functions (as `starts_with`) or  ina vector `c()`."))}


 #It is not common to have Condition column in Proline output. Account for its absence.

if(!missing(cond_col)){
  df <- df |>
    dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}, {{cond_col}})}else{
      df <- df |>
        dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}})
    }



  #if data in wide format, one must provide df_meta containing the same names as columns
  #with intensity values  as "Sample Names" with respective condition as "Condition"
  if(format == "wide"){
    if(!missing(cond_col)){

    df <-  df |>
    tidyr::pivot_longer(cols = -c({{seq_col}}, {{ptm_col}}),
                        names_to = "SampleName",
                        values_to = "intensity")
    }else{

      df <-  df |>
        tidyr::pivot_longer(cols = -c({{seq_col}}, {{ptm_col}}, {{cond_col}}),
                            names_to = "SampleName",
                            values_to = "intensity")

    }

    if(!missing(df_meta)) df <- dplyr::left_join(df, df_meta, by = "SampleName")

    }


  # This is to color the boxplots and points according to the condition and avoid error if no condition.
  if(!"Condition" %in% colnames(df) && missing(cond_col)) {COND_COL = FALSE
  cli::cli_alert_warning(cli::col_red("No column named 'Condition' was provided."))
  }else{COND_COL = TRUE}




  # split into two dataframes
    #replace PTM with 'xx'

  unique_seq <- df |> dplyr::select({{seq_col}}) |> unique() |>  dplyr::pull()

  plots <- list()

  for(seq in unique_seq){

  df_seq <- df |> dplyr::filter({{seq_col}} == seq)

  df_me1 <- df_seq |>
    dplyr::filter(stringr::str_detect({{ptm_col}}, {{me1}})) |>
    dplyr::mutate(PTMx = stringr::str_replace_all({{ptm_col}}, {{me1}}, "xx"), .after = {{ptm_col}})


  df_me1_label <- df_seq |>
    dplyr::filter(stringr::str_detect({{ptm_col}}, {{me1_label}})) |>
    dplyr::mutate(PTMx = stringr::str_replace_all({{ptm_col}}, {{me1_label}}, "xx"), .after = {{ptm_col}})

  #return(lst(df_me1, df_me1_label))

    #find common PTMs to compare

  me1_label_doublets <- dplyr::inner_join(unique(df_me1['PTMx']),
                                          unique(df_me1_label['PTMx']),
                                          by = "PTMx")

  #Join two dataframes in one dataframe
  df_me1 <- dplyr::right_join(df_me1,  me1_label_doublets, by = "PTMx")
  df_me1_label <- dplyr::right_join(df_me1_label,  me1_label_doublets, by = "PTMx")


  df_doublets <- dplyr::bind_rows(df_me1_label,df_me1)





  }



if(COND_COL == TRUE && !missing(cond_col)){

  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = {{cond_col}})) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank())



  plots[[seq]] <- p

}else if (COND_COL == TRUE && missing(cond_col)){
  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = Condition)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank())


  plots[[seq]] <- p
}else{

  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = "No condition")) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank())


  plots[[seq]] <- p
}




  if(save_plot){
   pdf("me1_labelling_efficiency_plots.pdf")
    for (plot_name in names(plots)) {
      cat(paste0('Saving plot: ', plot_name))
      print(plots[[plot_name]])
    }
    invisible(dev.off())

  }else{

    return(plots)}

}
