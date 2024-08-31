#' @title Intensity vs PTM grouped bar plot with jittered points
#' @description
#' Individual measurements will be plotted as data points. Their median (or mean) will be represented as a bar. Each bar will be colored according to the condition.
#'
#' @param dataset A dataframe in long format with at least 4 columns: Intensity, PTM, sequence (or sequence label) and Condition. check the arguemnts below.
#' A 5th column can be provided to be the plot title.
#' @param x_axis x variable (PTM column)
#' @param y_axis y variable (intensity column). If already values are percentage, set \code{scale} to 1.
#' @param condition The condition column (WT vs disease, concentration, ...)
#' @param id_col unique ID column such as sequence or sequence label
#' @param title_plot (\code{optional})A column with values to be the plot_title (optional).
#' @param fun median (\code{default}) or mean. This will be the height of the bar.
#' @param scale 100 (\code{default}). If you want to keep values as they are use 1. No other values are allowed.
#' @param cond_order (optional). A character vector containing conditions according to which bars will be ordered in the plot.
#' @param save_plot (\code{logical}; \code{optional})
#'
#' @import dplyr
#' @importFrom stats reorder
#' @importFrom scales label_percent
#' @importFrom rlang is_empty quo_name enquo
#' @importFrom stringr str_glue
#' @importFrom cli cli_abort cli_inform
#' @importFrom tidyr nest
#' @importFrom purrr map map2
#' @import ggplot2
#'
#' @return A dataframe grouped by `id_col` (and `plot_title` if passed) with a nested list-column harboring the generated bar plots
#' (representing either `mean` or `median`) with jitter points (corresponding to individual measurements per `condition`).
#' @export

plot_jitterbarIntvsPTM <- function(dataset,
                             x_axis, #PTM
                             y_axis, #Intensity
                             condition,  #variable on which comparison is done
                             id_col,  #could be sequence or sequence label
                             plot_title= NULL, #optional: sould be stripped sequence
                             fun = c("median", "mean"),
                             scale = c(100, 1),
                             cond_order= NULL,
                             save_plot = FALSE
){

# Check inputs------------
  fun <- match.arg(fun)
  stopifnot("Error: `scale` must be either 1 or 100." = scale %in% c(1, 100))
  if(missing(condition)){cli::cli_abort('`Condition` column is missing')}
  if(missing(id_col)){cli::cli_abort('`id_col` column is missing')}

  # prepare the dataset ------------

  dataset <- dataset |>
    dplyr::mutate({{condition}} := factor({{condition}}, levels = assort_cond(dataset,
                                                                                    condition_col = {{condition}},
                                                                                    cond_order = cond_order))) |>
    tidyr::nest(data = -c({{id_col}}, {{plot_title}})) |>
    dplyr::mutate(data = purrr::map(data, ~ .x %>%
                               dplyr::mutate(size = dplyr::n_distinct(.x[[rlang::quo_name(rlang::enquo(x_axis))]]))


                               )) |>
    dplyr::mutate(plots = purrr::walk2(.x= {{id_col}},
                        .y= data,
                        .f = ~ plotjit(id_col= .x,
                                       dataset = .y,
                                       x_axis = {{x_axis}},
                                       y_axis = {{y_axis}},
                                       condition = {{condition}},
                                       plot_title = {{plot_title}},
                                       fun= fun,
                                       scale= scale,
                                       save_plot = save_plot)))

}


#' @noRd

plotjit <- function(dataset,
                      x_axis,
                      y_axis, #Intensity
                      condition,  #variable on which comparison is done
                      id_col,  #could be sequence or sequence label
                      plot_title= NULL, #optional: sould be stripped sequence
                      fun = c("median", "mean"),
                      scale,
                      save_plot){




p <- ggplot2::ggplot(dataset, ggplot2::aes(x = stats::reorder({{x_axis}},{{y_axis}}),
                                           y = {{y_axis}},
                                         color = {{condition}},
                                          fill = {{condition}})) +

    #specify if mean or median is to be calacualted. also geom could be changed.
    ggplot2::stat_summary(fun = fun, show.legend =  TRUE, geom = "bar",
                          color = "black", width = 0.5,
                          position = ggplot2::position_dodge2(preserve = "single"),
                          alpha = 0.5) +
    ggplot2::geom_jitter(show.legend = FALSE,
                         position = ggplot2::position_jitterdodge(dodge.width = 0.5,
                         jitter.width = 0.1),
                         shape= 21, color = "black", size = ifelse(unique(dataset[['size']]> 10),
                                                                          unique(dataset[['size']]/10),
                                                                                 3), alpha = 1 ) +

    ggplot2::scale_y_continuous(labels = scales::label_percent(scale = scale)) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set1") +

    ggplot2::labs(x= "",
                  y= stringr::str_glue('% of all variably modified forms of {id_col}'),
                  title = ifelse(is.null(plot_title), "", {{plot_title}})) +

  ggplot2::guides(fill  = ggplot2::guide_legend(position = "inside")) +

  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x  = ggplot2::element_blank(),
    panel.grid.major.y  = ggplot2::element_line(color = "grey70"),
    axis.text.y = ggplot2::element_text(
      size = 14,
      margin = ggplot2::margin(
        t = 0,
        r = 20,
        b = 0,
        l = 20
      )
    ),
    axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
    axis.text.x = ggplot2::element_text(
      size = 12,
      angle = 60,
      colour = "black",
      hjust = 0.7
    ),
    axis.line.x = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0.5, 0.5, 0.2, 0.5), "cm") ,
    legend.text = ggplot2::element_text(face = "bold", size = 13),
    legend.position.inside = c(0.01, 0.95),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    title = ggplot2::element_text(size = 16)
  )


   cli::cli_inform("Plotting: {id_col}")



 if(save_plot){
     ggplot2::ggsave(filename = stringr::str_glue("{id_col}.png"),
            dpi = 300,
            height = 7,
            width = 10,
            units =  "in",
            bg = "white")
   }

   return(p)

}


#' @noRd



assort_cond <- function(data, condition_col,  cond_order) {

  conditions <- data |> dplyr::select({{condition_col}}) |> dplyr::pull()

  if (!is.null(cond_order)) {
    check_diff = base::setdiff(cond_order,unique(conditions)) #the same function found in tidyverse (lubridate or dplyr)

    if(!rlang::is_empty(check_diff)){
      cli::cli_abort(c( 'x' = 'The provided condition in `cond_order` does not match the ones in your dataset.',
                        'i' = 'Check the following condition(s): {check_diff}'))}else{
    sorted_levels <-  cond_order
                                               }
  } else {

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


