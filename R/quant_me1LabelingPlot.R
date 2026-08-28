
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
#' @param me1 Monomethyl representation inside the `ptm_col`column (i.e. default: "(?<![A-Za-z_])me1(?![A-Za-z])" or 'Methyl', ...)
#' @param me1_label Monomethyl+label representation inside the `ptm_col`column (i.e. default: 'bu' or 'Butyryl, ....)
#' @param isNormalized bool; FALSE (default).
#' @param scale 1 or 100 (default). Multiply (normalized) intensity values by a scalar.
#' @param format if set to `wide` (i.e. multiple samples), data will be reshaped into long format before plotting. Check `df_meta`. If `long`, the `int_col` will be renamed `intensity`.
#' @param df_meta (`optional`). If df is in `wide` format, then one can provide df_meta containing at least two columns (`SampleName` and `Condition`). This will be used during plotting.
#' `SampleNames` MUST match EXACTLY the names of the `int_col`columns.
#' @param save_plot (`logical`). If `TRUE`, all the plots will be saved in one pdf. If `FALSE` (`default`), the plot(s) will be shown within the Plots pane.)
#' @param save_data (`logical`). If `TRUE`, data containing all me1 or me1-labeled instances are saved in a csv file.
#' @param plotfile_name (`character`). Name of the pdf file containing all the comparison box plots.
#' @param datafile_name (`character`). Name of the csv file containing data of the plotted boxplots.
#' @param output_dir (working directory (default)). A folder name to output the plot df file and the csv data file.
#'
#' @importFrom dplyr select left_join filter mutate inner_join right_join pull bind_rows
#' @importFrom stringr str_detect str_replace_all str_wrap
#' @importFrom cli cli_alert_warning col_red cli_abort cli_alert_info
#' @importFrom rlang as_name ensym
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom readr write_csv
#' @importFrom grDevices dev.off pdf
#' @return A list of two: `plots` list containing boxplot with jittered points per sequence and `data` containing a tibble storing all the data used to generate these plots.
#' @export
#'

quant_me1LabelingPlot <- function(df,
                                  seq_col,
                                  ptm_col,
                                  int_col,
                                  cond_col= NULL,
                                  me1 = "(?<![A-Za-z_])me1(?![A-Za-z])",
                                  me1_label = "bu",
                                  isNormalized = FALSE,
                                  scale = 100,
                                  format= c('wide', 'long'),
                                  df_meta,
                                  save_plot = FALSE,
                                  save_data = FALSE,
                                  plotfile_name = 'plots_me1.pdf',
                                  datafile_name = 'data_me1.csv',
                                  output_dir= "."){
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

 df_check <- df |> dplyr::select(dplyr::all_of({{int_col}}))

 if (ncol(df_check) == 0) {
   cli::cli_abort(c('x' = "The provided column names in `int_col` argument are not found in `df`.",
                    'i' = "provide correct names using `tidyselect` functions (as `starts_with`) or  in a vector `c()`."))
   }


 if (format == "long") {
   df <- if (!is.null(cond_col)) {
     df |> dplyr::rename(intensity = {{int_col}}) |>  dplyr::select({{seq_col}}, {{ptm_col}}, {{cond_col}}, intensity)
   } else {
     df |> dplyr::rename(intensity = {{int_col}}) |> dplyr::select({{seq_col}}, {{ptm_col}}, intensity)
   }
 } else if (format == "wide") {
   df <- df |>
     dplyr::select({{seq_col}}, {{ptm_col}}, dplyr::all_of({{int_col}})) |>
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

if(isNormalized == FALSE){
 df <- quant_relIntensity(df = df,
                          select_cols = 'intensity',
                          grouping_var = dplyr::all_of(c('SampleName',
                                                         rlang::as_name(rlang::ensym(seq_col))))
 )}

 df <- df |> dplyr::mutate(intensity = scale*intensity)



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
data_plots <- list()

#loop over all unique sequences
 for(sq in unique_seq){

  df_seq <- df |> dplyr::filter({{seq_col}} == sq)


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

  df_me1 <- flexy_join(df_me1,  me1_label_doublets, by_col= "PTMx")
  df_me1_label <- flexy_join(df_me1_label,  me1_label_doublets, by_col= "PTMx")

  #Join me1 and me1-labed in one dataframe
  df_doublets <- dplyr::bind_rows(df_me1_label,df_me1)

# discard sequences with no me1 at all
if (nrow(df_doublets) > 0) {

  # color mapping depends on the COND_COL / cond_col state determined earlier
  # (see the "Condition" column detection block above)
  color_aes <- if (COND_COL == TRUE && !is.null(cond_col)) {
    rlang::quo({{ cond_col }})          # user explicitly passed a column
  } else if (COND_COL == TRUE && is.null(cond_col)) {
    rlang::quo(Condition)               # merged/found "Condition" column
  } else {
    rlang::quo("None")                  # no condition available, group as "None"
  }

   df_plot <- df_doublets |>
    tidyr::drop_na(intensity)

  p <- ggplot2::ggplot(data = df_plot, ggplot2::aes(y = intensity, x = {{ ptm_col }}, color = !!color_aes)) +
    ggplot2::geom_boxplot(
      width = 0.3,
      position = ggplot2::position_dodge(width = 0.4),
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(position = ggplot2::position_dodge2(width = 0.2)) +
    ggplot2::labs(title = sq, y = "% Relative Intensity") +
    ggplot2::scale_x_discrete(labels = wrap_labels) +
    ggplot2::facet_grid(~ PTMx, scales = "free_x", space = "free_x") +
    ggplot2::scale_color_viridis_d(option = "H") +
    ggplot2::theme_minimal(base_size = 12)+
    ggplot2::theme(
     #axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
     panel.spacing = ggplot2::unit(0.5, "lines")   # default is ~1 line;
    )

  plots[[sq]] <- p
  data_plots[[sq]] <- df_plot

}
 }


if (length(plots) == 0) {
  cli::cli_alert_warning("No plots were generated — check your input data.")
  return(NULL)
}

combined_data <- dplyr::bind_rows(data_plots, .id = "sequence")

if(save_plot|save_data){
output_dir <- if (is.null(output_dir)) getwd() else output_dir
if (!is.null(output_dir) && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cli::cli_inform('The folder {.path {output_dir}} is created.')
}
}

  if(save_plot){

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    if (!grepl("\\.pdf$", plotfile_name, ignore.case = TRUE)) {
      plotfile_name <- paste0(plotfile_name, ".pdf")
    }

    plot_path <- file.path(output_dir, plotfile_name)

    pdf(plot_path, width = 11, height = 9, author = 'histonePTM')
    on.exit(dev.off(), add = TRUE)
    for (plot_name in names(plots)) {
      cat(paste0('Saving plot: ', plot_name, "\n"))

      p <- plots[[plot_name]]

      # count facet panels to decide how cramped this plot is
      n_panels <- length(ggplot2::ggplot_build(p)$layout$panel_params)

      # scale text size down as panel count grows
      label_size <- dplyr::case_when(
        n_panels <= 2 ~ 9,
        n_panels <= 4 ~ 7.5,
        n_panels <= 6 ~ 6.5,
        TRUE          ~ 5.5
      )

      p <- p +
        ggplot2::theme(
          strip.text = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = label_size),
          plot.margin = ggplot2::margin(t = 5, r = 5, b = 40, l = 5)
        )

      print(p)
    }

    cli::cli_alert_success("Plots saved to {.path {plot_path}}")


  }

if (save_data) {

  if (!grepl("\\.csv$", datafile_name, ignore.case = TRUE)) {
    datafile_name <- paste0(datafile_name, ".csv")
  }
  data_path <- file.path(output_dir, datafile_name)
  readr::write_csv(combined_data, data_path)
  cli::cli_alert_success("Data saved to {.path {data_path}}")
}

return(invisible(list(plots = plots, data = combined_data)))
}



#' @noRd
wrap_by_hyphen <- function(label, max_width = 25) {
  # split, keeping the hyphen attached to the end of each chunk
  parts <- stringr::str_split(label, "(?<=[-\u2212])")[[1]]

  lines <- character(0)
  current <- ""
  for (p in parts) {
    candidate <- paste0(current, p)
    if (nchar(candidate) > max_width && current != "") {
      lines <- c(lines, current)
      current <- p
    } else {
      current <- candidate
    }
  }
  lines <- c(lines, current)
  paste(lines, collapse = "\n")
}

#' @noRd
wrap_labels <- function(labels, max_width = 17) {
  sapply(labels, function(label) {
    if (stringr::str_detect(label, "; ")) {                 # Proline-style semicolon lists
      stringr::str_replace_all(label, "; ", "\n")
    } else if (stringr::str_detect(label, "[-\u2212]")) {    # hyphen-chained labels (both - and −)
      wrap_by_hyphen(label, max_width)
    } else {                                                  # normal space-separated text
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
