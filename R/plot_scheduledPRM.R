#' plot_scheduledPRM
#'
#' @param df A \code{data.frame} or \code{tibble} having at least two columns: mz and rt.
#' @description
#' A function to plot scheduled PRM graph to check how many precursor ions are concurrently
#' monitored over the same time by varying the \code{rt_window} size.
#'
#' @param rt_col retention time column (min)
#' @param mz_col mass-over-charge column
#' @param rt_window The retention time window in min (a 5-min window means +-2.5 min).
#' multiple values can be submitted as a vector (e.g \code{c(5, 7)}).
#' @param peak_width The width of the elution peak in min (Default is 0.3).
#' @param save_plot a \code{logical} value to either save \code{TRUE} or to only view the plot \code{FALSE}.
#'  default value is \code{TRUE}.
#' @import rlang
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @import grid
#' @import tidyr
#' @import cli
#'
#' @examples
#' plot_scheduledPRM(mtcars, mpg, cyl, c(5, 10), save_plot = TRUE)
#'
#' @return plot.
#' @export
plot_scheduledPRM <- function(df, rt_col, mz_col, rt_window, peak_width =0.3, save_plot = FALSE) {

  df <- tibble::as_tibble(df)

  elution_peak_width_min = peak_width

#find min and max
  range <- df |> dplyr::summarise(min = min({{rt_col}}, na.rm = TRUE),
                                  max = max({{rt_col}}, na.rm = TRUE))

  tr_min <- range |>  purrr::pluck(min)
  tr_max <-  range |>  purrr::pluck(max)

  #set the bins according to the defined elution peak width
  bins <- seq(from = tr_min, to= tr_max, by = elution_peak_width_min)

  if(bins[length(bins)] < tr_max){
    bins <- append(bins, tr_max)}

  df1 <- rt_window |> (\(rt) {
    purrr::map_dfr(.x = rt,  ~{df |>  dplyr::select({{rt_col}}, {{mz_col}}) |>

        dplyr::mutate(tr_start = {{rt_col}} - .x/2,
                      tr_end   = {{rt_col}} + .x/2,
                      tr_win = .x
        )})
  })()


  #make all combination between generated bins
  df2 <- tidyr::crossing(df1, bins)

  #if precursor mass lies within certain bin add 1, otherwise add 0

  df3 <- df2 |>
    dplyr::mutate(see_if = dplyr::if_else(dplyr::between(bins, tr_start, tr_end),
                                          1,
                                          0), .by = tr_win)  |>
    tidyr::pivot_wider(names_from = bins, values_from = see_if) |>
    dplyr::select(-c(mpg, cyl, tr_start, tr_end)) |>
    dplyr::group_by(tr_win) |>
    dplyr::summarise_all(sum) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(
      cols = -tr_win,
      names_to = "time",
      values_to = "num_prec",
      names_transform = as.double,
      values_transform = as.integer
    ) |>
    dplyr::mutate(tr_win = paste0(tr_win, " min"))

  #draw the plot

  df3 |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = num_prec)) +
    ggplot2::geom_line(linewidth = 1.1,  ggplot2::aes(color = tr_win)) +
    ggplot2::labs(x = "Scheduled Time (min)", y = "Concurrent Precursors", color = "") +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linewidth = 2))) +
    ggplot2::coord_cartesian(clip = 'off') +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linetype = "dashed"),
      legend.key.width = grid::unit(1.5, "cm"),
      panel.grid.major.x =  ggplot2::element_blank(),
      legend.direction = "horizontal",
      legend.position = "top"
    )

  if(isTRUE(save_plot)){

    ggplot2::ggsave(stringr::str_glue("prm_scheduled.png"), height = 7, width =10, dpi = 300)
    cli::cli_alert_success(text = "Scheduled PRM plot is saved successfully.")
  }

}



