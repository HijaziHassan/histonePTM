#' @title Intensity vs PTM grouped bar plot with jittered points
#' @description
#' Individual measurements are plotted as points with a bar representing
#' the median (or mean). Bars are grouped/colored per condition.
#'
#' @param dataset A dataframe in long format with at least 4 columns: Intensity, PTM, sequence (or sequence label. see `id_col`) and Condition.
#' A 5th column can be provided to be the plot title.
#' @param x_axis x_axis variable (PTM column)
#' @param y_axis y_axis variable (intensity column). If already values are percentage, set \code{scale} to 1.
#' @param condition The condition column (WT vs disease, concentration, ...)
#' @param id_col Unique ID column such as sequence or sequence label
#' @param plot_title (optional) A column with values to be the plot_title.
#' @param max_cutoff (optional) Maximum value above which values will be filtered out. It is added to draw separately low abundant PTMs.
#' @param fun "mean" (default) or "median" - the height of the bar.
#' @param error_type one of  "none" (default), "CI" (confidence interval), "SE" (standard error), "SD" (standard deviation).
#' @param conf_level Confidence level (default 0.95)
#' @param scale 100 (default). If you want to keep values as they are, use 1. No other values are allowed.
#' @param cond_order (optional). Character vector of condition order.
#' @param save_plot Logical. If `TRUE`, saves the plot to disk. Default is `FALSE`.
#' @param output_dir (character; Optional) The output directory for the plots to be saved in. Default to working directory if unassigned.
#' @param base_point_size Base point size (will be scaled down for many categories)
#' @param base_bar_width Base bar width (will be scaled down for many categories)
#' @param base_text_size Base text size (will be scaled down for many categories)
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
                             output_dir = ".",
                             base_point_size = 2.5,
                             base_bar_width = 0.7,
                             base_text_size = 13
){


# Check inputs------------
  fun <- match.arg(fun)
  error_type <- match.arg(error_type)
  # output_dir <- if (is.null(output_dir)) getwd() else output_dir
  # if (!is.null(output_dir) && !dir.exists(output_dir)) {
  #   dir.create(output_dir, recursive = TRUE)
  # }


  stopifnot("Error: `scale` must be either 1 or 100." = scale %in% c(1, 100))
  if(missing(condition)){cli::cli_abort('`Condition` column is missing')}
  if(missing(id_col)){cli::cli_abort('`id_col` column is missing')}
  #If PTM column used in grouping, it will no more be available in the nested dataframe
      # ==> Check if x_axis and plot_title arguments are not identical

  # Ensures x_axis and y_axis are symbols (i.e. column names)
  x_axis <- rlang::ensym(x_axis)
  y_axis <- rlang::ensym(y_axis)

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



  # Apply max_cutoff filter
  if (is.numeric(max_cutoff)) {
    dataset <- dataset |> dplyr::filter({{ y_axis }} <= max_cutoff)
  } else if (!is.null(max_cutoff)) {
    cli::cli_alert_danger("`max_cutoff` must be numeric – ignoring.")
  }



  # Handle missing condition values
    condition_col <- dplyr::pull(dataset, {{condition}})

    if (all(is.na(condition_col))) {
      dataset <- dataset |> dplyr::mutate(!!rlang::ensym(condition) := "cond_1")
      cli::cli_alert_warning("All values in the column were NA. Replaced with 'cond_1'.")
    } else if (any(is.na(condition_col))) {
      dataset <- dataset |> dplyr::mutate(!!rlang::ensym(condition) := tidyr::replace_na(condition_col, "cond_2"))
      cli::cli_alert_warning("Some values in the column were NA - Replaced with 'cond_2'.")
    }



  # Factor conditions and clean x_axis
  dataset <- dataset |>
    dplyr::mutate({{condition}} := factor({{condition}},
                                          levels = assort_cond(dataset,
                                                               condition_col = {{condition}},
                                                               cond_order = cond_order)),
                  {{ x_axis }} := tidyr::replace_na(as.character({{ x_axis }}), "no_PTM")
                  )

  global_cond_levels <- levels(dataset[[rlang::as_label(rlang::enquo(condition))]])


  id_col <- rlang::enquo(id_col)
  plot_title <- rlang::enquo(plot_title)

  # ---- Nest and plot ----

  plot_title_label <- rlang::as_label(plot_title)


  # Check if plot_title is NULL or provided as a string (quoted)
  if (rlang::quo_is_null(plot_title) || !plot_title_label %in% names(dataset)) {
    dataset <- dataset |>
      tidyr::nest(data = -!!id_col) |>
      dplyr::mutate(
        n_categories = purrr::map_int(data, ~ dplyr::n_distinct(.x[[rlang::quo_name(x_axis)]])),
        plots = purrr::pmap(
          list(!!id_col, data, n_categories),
          function(id, dat, ncat) {
            plotjit(
              id_col = id,
              dataset = dat,
              x_axis = {{ x_axis }},
              y_axis = {{ y_axis }},
              condition = {{ condition }},
              plot_title = NULL,
              fun = fun,
              error_type = error_type,
              conf_level = conf_level,
              scale = scale,
              global_cond_levels = global_cond_levels,
              n_categories = ncat,
              base_point_size = base_point_size,
              base_bar_width = base_bar_width,
              base_text_size = base_text_size,
              save_plot = save_plot,
              output_dir = output_dir
            )
          }
        )
      )
  } else {
      # With plot_title – nest by id_col and plot_title
      dataset <- dataset |>
        tidyr::nest(data = -c(!!id_col, !!plot_title)) |>
        dplyr::mutate(
          n_categories = purrr::map_int(data, ~ dplyr::n_distinct(.x[[rlang::quo_name(x_axis)]])),
          plots = purrr::pmap(
            list(!!id_col, data, !!plot_title, n_categories),
            function(id, dat, title, ncat) {
              plotjit(
                id_col = id,
                dataset = dat,
                x_axis = {{ x_axis }},
                y_axis = {{ y_axis }},
                condition = {{ condition }},
                plot_title = title,
                fun = fun,
                error_type = error_type,
                conf_level = conf_level,
                scale = scale,
                global_cond_levels = global_cond_levels,
                n_categories = ncat,
                base_point_size = base_point_size,
                base_bar_width = base_bar_width,
                base_text_size = base_text_size,
                save_plot = save_plot,
                output_dir = output_dir
              )
            }
          )
        )
    }

  return(dataset)
}

#' @noRd
plotjit <- function(dataset,
                    x_axis,
                    y_axis,
                    condition,
                    id_col,
                    plot_title = NULL,
                    max_cutoff = NULL,
                    fun = c("mean", "median"),
                    error_type = c("none", "CI", "SE", "SD"),
                    conf_level = 0.95,
                    scale = 100,
                    save_plot = FALSE,
                    output_dir = '.',
                    global_cond_levels = NULL,
                    n_categories = NULL,
                    base_point_size = 2.5,
                    base_bar_width = 0.7,
                    base_text_size = 13) {

  fun <- match.arg(fun)
  error_type <- match.arg(error_type)

  # ---- Initial empty dataset guard ----
  if (nrow(dataset) == 0) {
    cli::cli_alert_warning("No data for {id_col} – skipping plot.")
    return(ggplot2::ggplot() +
             ggplot2::theme_void() +
             ggplot2::labs(title = "No data"))
  }

  # ---- Compute n_categories if not provided ----
  if (is.null(n_categories)) {
    n_categories <- dataset |>
      dplyr::pull({{ x_axis }}) |>
      unique() |>
      length()
  }

  # ---- Adaptive sizing ----
  scale_factor <- n_categories^(-0.15)

  # ---- Clean data (filter out NAs and groups with < 1) ----
  dataset <- dataset |>
    tidyr::drop_na({{ y_axis }}) |>
    dplyr::group_by({{ x_axis }}, {{ condition }}) |>
    dplyr::filter(dplyr::n() >= 1) |>
    dplyr::ungroup()

  # ---- Re-check after filtering ----
  if (nrow(dataset) == 0) {
    cli::cli_alert_warning("No data remaining for {id_col} after filtering – skipping plot.")
    return(ggplot2::ggplot() +
             ggplot2::theme_void() +
             ggplot2::labs(title = "No data after filtering"))
  }

  # Clamp scale_factor to prevent extreme values
  scale_factor <- max(0.95, min(1.2, scale_factor))

  point_size <- base_point_size * scale_factor
  bar_width <- base_bar_width * scale_factor
  x_axis_text_size <- base_text_size * scale_factor
  x_axis_text_size <- max(x_axis_text_size, 8)   # never smaller than 8 pt

  # Angle for x-axis labels: 90° if many categories with long labels
  max_label_len <- dataset |>
    dplyr::pull({{ x_axis }}) |>
    as.character() |>
    nchar() |>
    max(na.rm = TRUE)   # <-- added na.rm = TRUE for safety

  angle <- if (n_categories > 10 && max_label_len > 11) 90 else 60

  # ---- Compute summary statistics for error bars if applicable ----
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

  # Build named color/fill vectors from the global palette
  n_cond <- length(global_cond_levels)
  palette_colors <- if (!is.null(global_cond_levels) && n_cond > 0) {
    setNames(RColorBrewer::brewer.pal(max(3, n_cond), "Set1")[seq_len(n_cond)],
             global_cond_levels)
  } else NULL

  # Create the plot
  p <- ggplot2::ggplot(dataset, ggplot2::aes(
    x = stats::reorder({{ x_axis }}, {{ y_axis }}),
    y = {{ y_axis }},
    fill = {{ condition }}
  )) +
    ggplot2::stat_summary(
      fun = fun,
      show.legend = TRUE,
      geom = "bar",
      width = bar_width,
      position = ggplot2::position_dodge2(preserve = "single"),
      alpha = 0.5,
      aes(fill = {{ condition }})
    ) +
    ggplot2::geom_jitter(
      aes(color = {{ condition }}),
      position = ggplot2::position_jitterdodge(
        dodge.width = bar_width,
        jitter.width = 0.1
      ),
      shape = 19,
      size = point_size,
      alpha = 1,
      show.legend = FALSE
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
      width = 0.0,
      size = 0.9,
      position = ggplot2::position_dodge(width = bar_width),
      inherit.aes = FALSE,
      color = "black"
    )
  }

  # Customize scales & labels
  p <- p +
    ggplot2::scale_y_continuous(labels = scales::label_percent(scale = scale)) +
    ggplot2::scale_colour_manual(values = palette_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = palette_colors, drop = FALSE) +
    ggplot2::labs(
      x = "",
      y = "% of variably modified forms",
      title = if (is.null(plot_title)) "" else plot_title
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey70"),
      axis.text.y = ggplot2::element_text(
        size = base_text_size,
        face = "bold",
        margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
      ),
      axis.text.x = ggplot2::element_text(
        size = x_axis_text_size,
        face = "bold",
        angle = angle,
        colour = "black",
        hjust = if (angle == 90) 1 else 0.7,
        vjust = if (angle == 90) 0.5 else 0.85
      ),
      axis.title.y = ggplot2::element_text(
        size = base_text_size,
        face = "bold",
        margin = ggplot2::margin(r = 3)
      ),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.2, 0.5), "cm"),
      legend.text = ggplot2::element_text(face = "bold", size = base_text_size),
      legend.position = "inside",
      legend.position.inside = c(0.01, 0.95),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      title = ggplot2::element_text(size = base_text_size)
    )

  # ---- Save plots ----
  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    plot_width <- 6 + n_categories * 0.25
    plot_height <- 6 + n_categories * 0.1
    ggplot2::ggsave(
      filename = stringr::str_glue("{id_col}.png"),
      path = output_dir,
      dpi = 300,
      height = min(plot_height, 14),
      width = min(plot_width, 20),
      units = "in",
      bg = "white"
    )
  }

  return(p)
}



#' @noRd

assort_cond <- function(data, condition_col,  cond_order) {

  conditions <- data |> dplyr::pull({{condition_col}})

  if (!is.null(cond_order)) {
    check_diff = base::setdiff(cond_order, unique(conditions))

    if(!rlang::is_empty(check_diff)){
      cli::cli_abort(c( 'x' = 'The provided conditions in `cond_order` do not match the ones in your dataset.',
                        'i' = 'Check the following condition(s): {check_diff}'
                        ))
      } else {
        return(cond_order)
      }
  }

    #to order conditions based on the numeric values (so 12 mM does not come before 5 mM for e.g.)
    unique_levels <- unique(conditions)
    numeric_values <- as.numeric(gsub("[^0-9.]", "", unique_levels))
    if (any(!is.na(numeric_values))) {
      return(unique_levels[order(numeric_values)])
    } else {
      return(sort(unique_levels))
    }
  }
