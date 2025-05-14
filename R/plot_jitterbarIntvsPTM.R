#' @title Intensity vs PTM grouped bar plot with jittered points
#' @description
#' Individual measurements will be plotted as data points. Their median (or mean) will be represented as a bar. Each bar will be colored according to the condition.
#'
#' @param dataset A dataframe in long format with at least 4 columns: Intensity, PTM, sequence (or sequence label) and Condition. check the arguments below.
#' A 5th column can be provided to be the plot title.
#' @param x_axis x variable (PTM column)
#' @param y_axis y variable (intensity column). If already values are percentage, set \code{scale} to 1.
#' @param condition The condition column (WT vs disease, concentration, ...)
#' @param id_col unique ID column such as sequence or sequence label
#' @param plot_title (\code{optional})A column with values to be the plot_title (optional).
#' @param max_cutoff (\code{optional; numeric} A the maximum value above which values will be filtered out. It is added to draw separately low abundant PTMs.
#' @param fun mean (\code{default}) or median. This will be the height of the bar.
#' @param error_type one of "CI" (confidence interval), "SE" (standard error), "SD" (standard deviation).
#' @param conf_level Confidence level (e.g. 0.95 (default), 0.99, etc...)
#' @param scale 100 (\code{default}). If you want to keep values as they are use 1. No other values are allowed.
#' @param cond_order (optional). A character vector containing conditions according to which bars will be ordered in the plot.
#' @param save_plot (\code{logical}; FALSE \code{default})
#' @param output_dir (character; Optional) The output directory for the plots to be saved in. Default to working directory if unassigned.
#'
#' @importFrom dplyr select pull n_distinct mutate filter
#' @importFrom stats reorder median qt
#' @importFrom scales label_percent
#' @importFrom rlang is_empty quo_name enquo ensym quo_is_null
#' @importFrom stringr str_glue
#' @importFrom cli cli_abort cli_inform cli_alert_danger cli_alert_warning
#' @importFrom tidyr nest drop_na replace_na
#' @importFrom purrr map map2
#' @import ggplot2
#'
#'
#'
#' @return A dataframe grouped by `id_col` (and `plot_title` if passed) with a nested list-column harboring the generated bar plots
#' (representing either `mean` or `median`) with jitter points (corresponding to individual measurements per `condition`).
#' @export

plot_jitterbarIntvsPTM <- function(dataset,
                             x_axis, #PTM
                             y_axis, #Intensity
                             condition,  #variable on which comparison is done
                             id_col,  #could be sequence or sequence label
                             plot_title= NULL, #optional: e.g. stripped sequence
                             max_cutoff = NULL, #optional
                             fun = c("mean", "median"),
                             error_type = c("none", "CI", "SE", "SD"),
                             conf_level = 0.95,
                             scale = 100,
                             cond_order= NULL,
                             save_plot = FALSE,
                             output_dir = NULL
){


# Check inputs------------
  fun <- match.arg(fun)
  output_dir <- if (is.null(output_dir)) getwd() else output_dir
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Ensures x_axis and y_axis are symbols (i.e. column names)
  x_axis <- rlang::ensym(x_axis)
  y_axis <- rlang::ensym(y_axis)

  stopifnot("Error: `scale` must be either 1 or 100." = scale %in% c(1, 100))
  if(missing(condition)){cli::cli_abort('`Condition` column is missing')}
  if(missing(id_col)){cli::cli_abort('`id_col` column is missing')}
  #If PTM column used in grouping, it will no more be available in the nested dataframe
      # ==> Check if x_axis and plot_title arguments are not identical

  # Capture the unevaluated plot_title expression
  plottitle_expr <- substitute(plot_title)

  # check if plot_title is a symbol
  if (is.symbol(plottitle_expr)) {
    # PTM and title cannot have the same column name.
    if (identical(x_axis, plottitle_expr)) {
      cli::cli_abort(c(
        "x" = "{.arg x_axis} and {.arg plot_title} are the same.",
        "i" = "Choose a different column or a static string in {.arg plot_title} to entitle your plot(s)."
      ))
    }
  }



  # prepare the dataset ------------


# make sure max_cutoff is numeric.

  if(is.symbol(substitute(max_cutoff))){

    max_cutoff_str <- deparse(substitute(max_cutoff))
    cli::cli_alert_danger("The intensity filter `{max_cutoff_str}` is ignored. {.arg max_cutoff} must be numeric.")

    }else if (is.character(max_cutoff)){

          cli::cli_alert_danger("The intensity filter {.arg {max_cutoff}} is ignored. {.arg max_cutoff} must be numeric.")

  }else if (is.numeric(max_cutoff)) {
       dataset <- dataset |>
         dplyr::filter({{y_axis}} <= max_cutoff)
   }




    condition_col <- dplyr::pull(dataset, {{condition}})

    if (all(is.na(condition_col))) {
      dataset <- dataset |> dplyr::mutate(!!rlang::ensym(condition) := "cond_1")
      cli::cli_alert_warning("All values in the column were NA. Replaced with 'cond_1'.")
    } else if (any(is.na(condition_col))) {
      dataset <- dataset |> dplyr::mutate(!!rlang::ensym(condition) := tidyr::replace_na(condition_col, "cond_2"))
      cli::cli_alert_warning("Some values in the column were NA. Replaced with 'cond_2'.")
    }




  dataset <- dataset |>
    dplyr::mutate({{condition}} := factor({{condition}},
                                          levels = assort_cond(dataset,
                                                               condition_col = {{condition}},
                                                               cond_order = cond_order)))


  id_col <- rlang::enquo(id_col)
  plot_title <- rlang::enquo(plot_title)
  plot_title_label <- rlang::as_label(plot_title)


  # Check if plot_title is NULL or provided as a string (quoted)
 if (rlang::quo_is_null(plot_title) || !plot_title_label %in% names(dataset)) {

   dataset <- dataset |>
     tidyr::nest(data = -!!id_col) |>
     dplyr::mutate(
       data = purrr::map(
         data,
         ~ .x |>
           dplyr::mutate(
             size = dplyr::n_distinct(.x[[rlang::quo_name(rlang::enquo(x_axis))]])
           )
       )
     )



   dataset <- dataset |>
     dplyr::mutate(
       plots = purrr::map2(
         .x = {{id_col}},
         .y = data,
         .f = ~ plotjit(
           id_col = .x,
           dataset = .y,
           x_axis = {{x_axis}},
           y_axis = {{y_axis}},
           condition = {{condition}},
           plot_title = {{plot_title}},
           max_cutoff = max_cutoff,
           fun = fun,
           error_type = error_type,
           conf_level = conf_level,
           scale = scale,
           save_plot = save_plot,
           output_dir = output_dir
         )
       )
     )




  } else {


    dataset <- dataset |>
      tidyr::nest(data = -c(!!id_col, !!plot_title)) |>
      dplyr::mutate(data = purrr::map(data, ~ .x |> #size to adjust the size of the jitter points
                                        dplyr::mutate(size = dplyr::n_distinct(.x[[rlang::quo_name(rlang::enquo(x_axis))]]))

      ))

    cols <- list(
      dataset |> dplyr::pull({{id_col}}),
      dataset |> dplyr::pull(data),
      dataset |> dplyr::pull({{plot_title}})
    )

    dataset <- dataset |>
      dplyr::mutate(
        plots = purrr::pmap(
          cols,
          .f = ~ plotjit(
            id_col = ..1,
            dataset = ..2,
            x_axis = {{x_axis}},
            y_axis = {{y_axis}},
            condition = {{condition}},
            plot_title = ..3,
            fun = fun,
            error_type = error_type,
            conf_level = conf_level,
            scale = scale,
            save_plot = save_plot,
            output_dir= output_dir
          )
        )
      )

}

}


#' @noRd

plotjit <- function(dataset,
                    x_axis,
                    y_axis,  # Intensity
                    condition,  # variable on which comparison is done
                    id_col,  # could be sequence or sequence label
                    plot_title = NULL,  # optional: should be stripped sequence
                    max_cutoff= NULL,
                    fun = c("mean", "median"),
                    error_type = c("none", "CI", "SE", "SD"),
                    conf_level = 0.95,
                    scale = 100,
                    save_plot = FALSE,
                    output_dir= NULL) {


  fun <- match.arg(fun)
  error_type <- match.arg(error_type)


  dataset <- dataset |>
    tidyr::drop_na({{ y_axis }}) |>
    dplyr::group_by({{ x_axis }}, {{ condition }}) |>
    dplyr::filter(dplyr::n() >= 1) |>
    dplyr::ungroup()

  #adjust angle based on x-axis labels
  angle = ifelse("size" %in% colnames(dataset) &&
                   unique(dataset[["size"]]) > 10 &&
                   any(nchar(dataset |> dplyr::pull({{x_axis}})) > 11),
                 90,
                 60
  )

  # Compute summary statistics for error bars if applicable
  if (error_type != 'none') {


    summary_stats <- dataset |>
      dplyr::group_by({{ x_axis }}, {{ condition }}) |>
      dplyr::summarize(
        y_value = if (fun == "mean") mean({{ y_axis }}, na.rm = TRUE) else median({{ y_axis }}, na.rm = TRUE),
        sd = sd({{ y_axis }}, na.rm = TRUE),
        n = dplyr::n(),
        se = sd / sqrt(n),
        ci = se * qt((1 + conf_level) / 2, df = n - 1),
        ymin = dplyr::case_when(
          error_type == "SD" ~ y_value - sd,
          error_type == "SE" ~ y_value - se,
          error_type == "CI" ~ y_value - ci,
          TRUE ~ NA_real_
        ),
        ymax = dplyr::case_when(
          error_type == "SD" ~ y_value + sd,
          error_type == "SE" ~ y_value + se,
          error_type == "CI" ~ y_value + ci,
          TRUE ~ NA_real_
        ),
        .groups = "drop"
      )
  }

  # Create the plot
  p <- ggplot2::ggplot(dataset, ggplot2::aes(
    x = stats::reorder({{ x_axis }}, {{ y_axis }}),
    y = {{ y_axis }},
   # color = {{ condition }},
    fill = {{ condition }}
  )) +
    ggplot2::stat_summary(
      fun = fun, show.legend = TRUE, geom = "bar",
      , width = 0.5,
      position = ggplot2::position_dodge2(preserve = "single"),
      alpha = 0.5
    )

  # Add error bars if specified
  if (error_type != 'none') {
    p <- p + ggplot2::geom_errorbar(
      data = summary_stats,
      ggplot2::aes(
        x = {{ x_axis }},
        ymin = ymin,
        ymax = ymax,
        group = {{ condition }}
      ),
      width = 0.0,  # Thin line
      size = 0.9,
      position = ggplot2::position_dodge(width = 0.5),
      inherit.aes = FALSE,  # Avoid inheriting aes from the main dataset
      color = "black"
    )

  }

  # Add jittered points
  p <- p + ggplot2::geom_jitter(aes(color = {{ condition }}),
    show.legend = FALSE,
    position = ggplot2::position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.1
    ),
    shape = 19,
    #color = "black",
    size = ifelse("size" %in% colnames(dataset) && unique(dataset[["size"]]) > 10,
                  2.75 / (1 + 0.05 * (unique(dataset[["size"]]) - 10)), #decay fn to adapt point size to # of PTM combinations
                  2.75),
    alpha = 1
  ) +

    # Customize scales
    ggplot2::scale_y_continuous(labels = scales::label_percent(scale = scale)) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set1") +

    # Add labels and title
    ggplot2::labs(
      x = "",
      y =  wrap_label(stringr::str_glue("% of all variably modified forms of {id_col}")),
      title = ifelse(is.null(plot_title), "", plot_title)
    ) +

    # Apply theme settings
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.grid.major.y  = ggplot2::element_line(color = "grey70"),
      axis.text.y = ggplot2::element_text(
        size = 15,
        face = "bold",
        margin = ggplot2::margin(
          t = 0,
          r = 20,
          b = 0,
          l = 20
        )
      ),
      axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 13,
        face = "bold",
        angle = angle,
        colour = "black",
        hjust = ifelse(angle == 90, 1, 0.7),
        vjust = ifelse(angle == 90, 0.5, 0.85)
      ),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.2, 0.5), "cm") ,
      legend.text = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "inside",
      legend.position.inside = c(0.01, 0.95),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      title = ggplot2::element_text(size = 17)
    )




  if (save_plot) {

    cli::cli_inform("Plotting: {id_col}")
    ggplot2::ggsave(
      filename = stringr::str_glue("{id_col}.png"),
      path = output_dir,
      dpi = 300,
      height = 7,
      width = 10,
      units = "in",
      bg = "white"
    )
  }

  return(p)
}






#' @noRd



assort_cond <- function(data, condition_col,  cond_order) {

  conditions <- data |> dplyr::select({{condition_col}}) |> dplyr::pull()

  if (!is.null(cond_order)) {
    check_diff = base::setdiff(cond_order, unique(conditions))

    if(!rlang::is_empty(check_diff)){
      cli::cli_abort(c( 'x' = 'The provided conditions in `cond_order` do not match the ones in your dataset.',
                        'i' = 'Check the following condition(s): {check_diff}'))}else{
    sorted_levels <-  cond_order
                                               }
  } else {

    #to order conditions based on the numeric values (so 12 mM does not come before 5 mM for e.g.)
    unique_levels <- unique(conditions)

    numeric_values <- as.numeric(gsub("[^0-9.]", "", unique_levels))

    if (any(!is.na(numeric_values))) {

      sorted_levels <- unique_levels[order(numeric_values)]

    } else {

    sorted_levels <- sort(unique_levels)

    }
  }

return(sorted_levels)

}



#' @noRd
wrap_label <- function(label, width = 35) {
if (nchar(label) > 47) {
  stringr::str_wrap(label, width = width)
} else {
  label
}
}
