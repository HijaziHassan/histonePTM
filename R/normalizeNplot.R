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
#' @param plot_title A string or a column name.
#' @param fun A statistical test ('mean' or 'median'). This will be the height of the bars.
#' @param error_type (optional) Argument to add error bars. It takes either "CI", "SD" or "SE" (check `plot_jitterbarIntvsPTM()`).
#' @param conf_level (optional) The confidence level (0.95 (default)). Only useful if you want to add error bars (see `error_type`).
#' @param cond_order (optional) A character vector containing the conditions in the preferred order. If you want a specific order of the condition.
#' @param save_plot (\code{bool}) TRUE (default) to save the plot(s). Otherwise, set it to FALSE.
#'
#' @importFrom dplyr mutate select left_join join_by starts_with
#' @importFrom tidyr pivot_longer
#' @importFrom openxlsx write.xlsx
#' @importFrom cli cli_abort
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
                           plot_title= NULL,
                           fun = "mean",
                           error_type = NULL,
                           conf_level = 0.95,
                           cond_order = NULL,
                           save_plot= TRUE){

labeling = match.arg(labeling)

  intensity_columns <- df_corrected |>
    dplyr::select({{int_cols}}) |>
    colnames()

  if(!all(intensity_columns %in% colnames(df_corrected))) {
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
    dplyr::select(PTM_unlabeled, {{seq_col}}, {{seq_stretch_col}}, {{int_cols}}) |>
    tidyr::pivot_longer(
      cols = {{int_cols}},
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
                         fun= fun,
                         error_type = error_type,
                         conf_level = conf_level,
                         scale= 100,
                         cond_order = cond_order,
                         plot_title = plot_title,
                         save_plot = save_plot)



}






