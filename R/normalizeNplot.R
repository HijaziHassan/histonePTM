#' Normalize and Plot
#' @description
#' There are cases where the output of `analyzeHistone()` function is not satisfactory. This happens when the software integrates the same peak for two nearly coeluting peptidoforms
#' as K18acK23un and K18unK23un. A need to re-normalize and to re-plot the data is the reason this function exists.
#'
#'
#' @param df_corrected A dataframe containing the corrected raw intensity values along with the necessary column (see below)
#' @param df_meta A dataframe containing at least the 'Condition' and 'SampleName'. The latter MUST be the same the names of the `int_cols`.
#' @param ptm_col A column containing PTM. The function `misc_clearLabeling()` according to `labeling` will be applied to remove any labeling.
#' @param seq_col A column containing stripped sequences. Will be used for normalization. Could also be used as `plot_title`.
#' @param seq_stretch_col A column containing the label of the sequences (will be used as title of the y axis).
#' @param int_cols Intensity columns. `tidyselect` functions can be used (e.g. starts_with('abundance_')) or character vector (e.g `c('col1', 'col2')`).
#' @param save_file (\code{bool}) TRUE (default) to save the normalized file.
#' @param labeling Labeling method (e.g. 'PA', 'TMA', 'PA_PIC', or 'none' (default)). Check `misc_clearLabeling()`.
#' @param rules (optional) Rules to unlabel PTM strings. Check `misc_clearLabeling()`.
#' @param ... Additional arguments passed to `plot_jitterbarIntvsPTM` other than the column names such as `fun`, `error_type`, `plot_title`, `save_plot` etc ... Revise the documentation for more details.
#'
#' @importFrom dplyr mutate select left_join join_by starts_with any_of all_of
#' @importFrom tidyr pivot_longer
#' @importFrom openxlsx write.xlsx
#' @importFrom cli cli_abort
#' @importFrom rlang is_empty
#' @return Same output as `plot_jitterbarIntvsPTM()`. A dataframe grouped by 'id_col' (and 'plot_title' if passed) with a nested list-column harboring the generated
#' bar plots (representing either 'mean' or 'median') with jitter points (corresponding to individual measurements per 'condition').
#' @export

normalizeNplot <- function(df_corrected,
                           df_meta,
                           ptm_col,
                           seq_col,
                           seq_stretch_col,
                           int_cols,
                           save_file = TRUE,
                           labeling= c('none', 'PA', 'TMA', 'PIC_PA'),
                           rules= NULL,
                           ...
                           ){

labeling = match.arg(labeling)

#check of all intensity column names are found in the dataset
  intensity_columns <- df_corrected |>
    dplyr::select(dplyr::any_of({{int_cols}})) |>
    colnames()

  if(!all(intensity_columns %in% colnames(df_corrected)) | rlang::is_empty(intensity_columns)) {
    cli::cli_abort(
      "The provided intensity column names in {.arg int_cols} are not found in {.arg df_corrected}."
    )
  }



  #normalize data

  df_norm <- quant_relIntensity(df = df_corrected,
                                select_cols = {{int_cols}},
                                grouping_var = {{seq_col}}) |>
    dplyr::mutate(PTM_unlabeled = misc_clearLabeling({{ptm_col}}, residue= 'remove', labeling = labeling, rules= rules),
                  PTM_unlabeled = stringr::str_remove_all(PTM_unlabeled, "\\*"), .after = {{ptm_col}})

  if(save_file){
    file = "corr_normalized.xlsx"
    df_norm |> openxlsx::write.xlsx(file = file)
    cli::cli_alert_success("The corrected and normalized values are saved as {file}.")
  }

  #transform data into long format. add a clean PTM_unlabeled column for plotting
  df_long <- df_norm |>
    dplyr::select(PTM_unlabeled, {{seq_col}}, {{seq_stretch_col}}, dplyr::all_of({{int_cols}})) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of({{int_cols}}),
      names_to = "sample",
      values_to = "intensity",
      values_drop_na = TRUE
    )


  #add meta data to use it for plotting
  df_labeled <- dplyr::left_join(df_long, df_meta, by = dplyr::join_by(sample == SampleName))



  #plot

 p <-  plot_jitterbarIntvsPTM(df_labeled,
                         x_axis = PTM_unlabeled,
                         y_axis= intensity,
                         condition = Condition,
                         id_col = {{seq_stretch_col}},
                         ...

                         )



}






