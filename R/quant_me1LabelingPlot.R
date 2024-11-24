
#' Assessing Monomethyl Derivatization
#'
#' @description
#' Evaluate the efficiency of monomethylated lysine derivarization per sequencexPTM combination.
#'
#' @param df Dataframe containing proteomic analysis dataset
#' @param seq_col Sequence column
#' @param ptm_col PTM column
#' @param int_col Intensity column(s)
#' @param cond_col (`optional`) Condition column (e.g. concentration) used for legend. If absent will be replaced by "None".
#' @param me1 Monomethyl representation inside the `ptm_col`column (i.e. 'me1' (`default`) or 'Methyl', ...)
#' @param me1_label Monomethyl+label representation inside the `ptm_col`column (i.e. 'bu' (`default`) or 'Butyryl, ....)
#' @param format if set to `wide` (most often), data will be reshaped into long format before plotting. Check `df_meta`.
#' @param df_meta (`optional`). If df is in `wide` format, then one can provide df_meta containing at least two columns (`SampleName` and `Condition`).
#' `SampleNames` MUST match exactly the names of the `int_col`columns.
#' @param save_plot (`logical`). If `TRUE`, all the plots will be saved in one pdf. If `FALSE` (`default`), the plot(s) will be shown within the Plots pane.)
#'
#' @importFrom dplyr select left_join filter mutate inner_join right_join pull bind_rows
#' @importFrom stringr str_detect str_replace_all str_wrap
#' @importFrom cli cli_alert_warning col_red cli_abort cli_alert_info
#' @importFrom rlang as_name ensym
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom grDevices dev.off pdf
#' @return Boxplot with jittered points to assess me1 labeling when using anhydrides to label lysines.
#' @export
#'

quant_me1LabelingPlot <- function(df,
                                  seq_col,
                                  ptm_col,
                                  int_col,
                                  cond_col= NULL,
                                  me1 = "me1",
                                  me1_label = "bu",
                                  format= c('wide', 'long'),
                                  df_meta,
                                  save_plot = FALSE){
  # Check for missing arguments
  missing_args <- c(
    if (missing(df)) "df",
    if (missing(seq_col)) "seq_col",
    if (missing(ptm_col)) "ptm_col",
    if (missing(int_col)) "int_col"
  )

  # Abort if any required arguments are missing
  if (length(missing_args) > 0) {
    if(length(missing_args) == 1){
    cli::cli_abort(c(
      "The following required argument is missing:",
      "x" = "{.arg {missing_args}}"
    ))} else {
      cli::cli_abort(c(
        "The following required arguments are missing:",
        "x" = "{.arg {missing_args}}"
      ))
    }
  }



format = match.arg(format)

  #Avoid forgetting or using wrong pattern inside tidyselect function to select intensity columns

 df_check <- df |> dplyr::select({{int_col}})

 if (ncol(df_check) == 0) {
   cli::cli_abort(c('x' = "The provided column names in `int_col` argument are not found in `df`.",
                    'i' = "provide correct names using `tidyselect` functions (as `starts_with`) or  in a vector `c()`."))
   }


 if (format == "long") {
   df <- if (!is.null(cond_col)) {
     df |> dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}, {{cond_col}})
   } else {
     df |> dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}})
   }
 } else if (format == "wide") {
   df <- df |>
     dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}) |>
     tidyr::pivot_longer(
       cols = -c({{seq_col}}, {{ptm_col}}),
       names_to = "SampleName",
       values_to = "intensity",
       values_drop_na = TRUE
     )
 }

 # remove sequences that are shared peptween histones
 df <- df |> dplyr::distinct(.keep_all = TRUE)

 # Join metadata if provided
 if (!missing(df_meta)) {
   df <- df |> dplyr::left_join(df_meta, by = "SampleName")

 }


 df <- quant_relIntensity(df = df,
                          select_cols = intensity,
                          grouping_var = dplyr::all_of(c('SampleName',
                                                         rlang::as_name(rlang::ensym(seq_col))))
 )

 df <- df |> dplyr::mutate(intensity = 100*intensity)



  # This is to color the boxplots and points according to the condition and avoid error if no condition.
 if (!"Condition" %in% colnames(df) && is.null(cond_col)) {
   COND_COL <- FALSE
   cli::cli_alert_warning(cli::col_red("No 'Condition' column is provided."))
   cli::cli_alert_info("If there are two or more conditions, they will be grouped for each modified sequence by 'none'.")

 } else if ("Condition" %in% colnames(df) && is.null(cond_col)) {
   COND_COL <- TRUE
   cli::cli_alert_info("'Condition' column from {.arg df_meta} is successfully merged.")

 } else {
   COND_COL <- TRUE
 }




#   split into two dataframes
#   replace PTM with 'xx'

unique_seq <- df |> dplyr::select({{seq_col}}) |> unique() |>  dplyr::pull()

plots <- list()

#loop over all unique sequences
 for(seq in unique_seq){

  df_seq <- df |> dplyr::filter({{seq_col}} == seq)


  df_me1 <- df_seq |>
    dplyr::filter(stringr::str_detect({{ptm_col}}, {{me1}})) |>
    dplyr::mutate(PTMx = stringr::str_replace_all({{ptm_col}}, {{me1}}, "xx"), .after = {{ptm_col}})


  df_me1_label <- df_seq |>
    dplyr::filter(stringr::str_detect({{ptm_col}}, {{me1_label}})) |>
    dplyr::mutate(PTMx = stringr::str_replace_all({{ptm_col}}, {{me1_label}}, "xx"), .after = {{ptm_col}})



    #find common PTMs to compare

  me1_label_doublets <- dplyr::inner_join(unique(df_me1['PTMx']),
                                          unique(df_me1_label['PTMx']),
                                          by = "PTMx")



  # keep those case where either me1 or me1-labeled is present but not both as well

  df_me1 <- flexy_join(df_me1,  me1_label_doublets, by = "PTMx")
  df_me1_label <- flexy_join(df_me1_label,  me1_label_doublets, by = "PTMx")

  #Join me1 and me1-labed in one dataframe
  df_doublets <- dplyr::bind_rows(df_me1_label,df_me1)


# discard sequences with no me1 at all
if(nrow(df_doublets)>0){

if(COND_COL == TRUE && !is.null(cond_col)){

  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = {{cond_col}})) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq, y= "% Relative Intensity")+
    ggplot2::scale_x_discrete(labels = wrap_labels)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))



  plots[[seq]] <- p

}else if (COND_COL == TRUE && is.null(cond_col)){

  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = Condition)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq, y= "% Relative Intensity")+
    ggplot2::scale_x_discrete(labels = wrap_labels)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


  plots[[seq]] <- p

}else{

  p <- df_doublets |>
    tidyr::drop_na(intensity) |>
    ggplot2::ggplot(ggplot2::aes(y = intensity, x = {{ptm_col}}, color = "None")) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width  = 0.5)) +
    ggplot2::labs(title= seq, y= "% Relative Intensity")+
    ggplot2::scale_x_discrete(labels = wrap_labels)+
    ggplot2::facet_wrap( ~ PTMx, scales = "free") +
    ggplot2::scale_color_viridis_d(option = "H")+
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


  plots[[seq]] <- p
}


}
  }

  if(save_plot){
   pdf("me1_labelling_efficiency_plots.pdf")
    for (plot_name in names(plots)) {
      cat(paste0('Saving plot: ', plot_name, "\n"))
      print(plots[[plot_name]])
    }
    invisible(dev.off())

  }else{

    return(plots)}

}



#' @noRd

wrap_labels <- function(labels, max_width = 25) {
  sapply(labels, function(label) {
    if (stringr::str_detect(label, "; ")) { #adapted to Proline software
      stringr::str_replace_all(label, "; ", "\n")
    } else {
      stringr::str_wrap(label, width = max_width)
    }
  })
}

#' @noRd

flexy_join <- function(df, ref_df, by_col) {
  # Perform right join
  joined_df <- dplyr::right_join(df, ref_df, by = by_col)

  # Check if the resulting dataframe is empty
  if (nrow(joined_df) == 0) {
    joined_df <- dplyr::left_join(df, ref_df, by = by_col)
  }

  return(joined_df)
}
