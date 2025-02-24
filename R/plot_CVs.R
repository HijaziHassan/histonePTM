#' Viusalizaing Coefficients of Variation (CVs)
#' @description
#' Plot distribution of CVs as violin plot and count IDs that have CV of 20 % and 5 % or less.
#'
#' @param df A dataframe.
#' @param cond_col The condition column.
#' @param cv_col The column containing the coefficients of variations
#' @param scale  1 or 100. If CVs are in percentage, use 1. If not, set it to 100.
#' @param save_plot Save plot (bool; FALSE (default))
#' @param output_dir Folder to save the plot in. Defualts to working directory.
#'
#' @importFrom dplyr summarise mutate n
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom scales label_percent
#' @importFrom cli cli_abort
#' @import ggplot2
#'
#' @return list of two ggplot2. the violin_plot and bar_plot of CVs
#' @export
#'
#'
#'
plot_CVs <- \(df, cond_col, cv_col, scale = 100, save_plot = FALSE, output_dir= NULL){

  output_dir <- if (is.null(output_dir)) getwd() else output_dir
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Validate missing arguments
  if (missing(df)) { cli::cli_abort('{.arg df} is missing.') }
  if (missing(cond_col)) { cli::cli_abort('{.arg cond_col} is missing.') }
  if (missing(cv_col)) { cli::cli_abort('{.arg cv_col} is missing.') }

  if (!deparse(substitute(cond_col)) %in% names(df)) {
    cli::cli_abort("Column {.val {deparse(substitute(cond_col))}} not found in the dataframe.")
  }
  if (!deparse(substitute(cv_col)) %in% names(df)) {
    cli::cli_abort("Column {.val {deparse(substitute(cv_col))}} not found in the dataframe.")
  }

  # Compute means per condition
  mean_values <- df |>
    dplyr::summarise(mean_CV = mean({{cv_col}}, na.rm = TRUE),
                     .by = {{cond_col}})

  # Compute bar graph data
  bar_data <- df |>
    tidyr::drop_na({{cv_col}}) |>
    dplyr::summarise(
      `Total ID's` = dplyr::n(),
      `CV < 20%` = sum({{cv_col}} <= (20 * (1/scale)), na.rm = TRUE),
      `CV < 5%` = sum({{cv_col}} <= (5 * (1/scale)), na.rm = TRUE),
      .by = {{cond_col}}
    ) |>
    tidyr::pivot_longer(cols = c(`Total ID's`, `CV < 20%`, `CV < 5%`),
                        names_to = "Metric", values_to = "Count") |>
    dplyr::mutate(Metric = factor(Metric, levels = c("Total ID's", "CV < 20%", "CV < 5%")))  # Ensure correct order

  # Violin + Box Plot
  violin_plot <- ggplot2::ggplot(df, ggplot2::aes(x= {{cond_col}}, y= {{cv_col}}, color = {{cond_col}})) +
    ggplot2::geom_violin(trim = FALSE) +
    ggplot2::geom_boxplot(linewidth = .8, width = .1, ggplot2::aes(color = {{cond_col}})) +
    ggplot2::geom_point(data = mean_values, aes(x = {{cond_col}}, y = mean_CV),
                        shape = 4, size = 3, color = "black") +  # Mean as X marker
    ggplot2::geom_text(data = mean_values, ggplot2::aes(x = {{cond_col}},
                                               y = mean_CV + 0.05,
                                               label = round(mean_CV, 2)),
                       color = "black", size = 5) +  # Mean text above box plot
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(y= "% CV", x= "") +
    ggplot2::scale_y_continuous(expand = c(0,0),
                                labels = scales::label_percent(scale = scale)) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 16, colour = "black"),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 1, colour = "grey90")
    )

  # Bar Plot (Bars closer, X-axis starts from 0)
  bar_plot <- ggplot2::ggplot(bar_data, ggplot2::aes(x = {{cond_col}}, y = Count, fill = Metric)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.32),
                      width = 0.3) +  # Closer bars
    ggplot2::geom_text(ggplot2::aes(label = Count),
                       position = ggplot2::position_dodge(width = 0.32),
                       vjust = -0.5, fontface = "bold", color = "black", size = 5) +  # Bold black labels on top
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(y = "# of Peptides", x = "") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_discrete(expand = c(0.3, 0))+
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(bar_data$Count) * 1.1)) +
    ggplot2::theme_classic(base_size = 16) +

    ggplot2::theme(
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 14, colour = "black", face = 'bold'),
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 16, colour = "black", face = 'bold'),
      axis.text.x = ggplot2::element_text(size = 16, colour = "black", face = 'bold'),
      axis.line.x =  ggplot2::element_line(linewidth = 1.2),
      panel.grid.major.y = ggplot2::element_line(linewidth = 1, colour = "grey90"),
      legend.margin = margin(t = -10, unit = "pt")
    )


  if (save_plot) {

    ggplot2::ggsave(
      violin_plot,
      filename = 'CV_Violinplot.png',
      path = output_dir,
      dpi = 300,
      height = 7,
      width = 6,
      units = "in",
      bg = "white"
    )

    ggplot2::ggsave(
      bar_plot,
      filename = 'CV_Barplot.png',
      path = output_dir,
      dpi = 300,
      height = 7,
      width = 6,
      units = "in",
      bg = "white"
    )


  }
  # Return both plots
  list(violin_plot = violin_plot, bar_plot = bar_plot)


}


