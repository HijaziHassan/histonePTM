


#' A ridgeline visualization of intensities per sample
#' @description
#' A plot the gives an eye-bird view of the intensity distribution in each sample.
#'
#' @param df Dataset in wide format.
#' @param int_col The intensity columns. Either a vector `c('column1', 'column2', ...)` or using tidyselect functions (e.g. starts_With('abundance'))
#' @param seq_col sequence column.
#' @param ptm_col modification column.
#' @param format 'wide' or 'long'.
#' @param save_plot \code{bool}, FALSE (default).
#' @importFrom ggridges geom_density_ridges_gradient
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select distinct mutate everything any_of
#'
#' @return ridgeline plot
#' @export
#'


plot_intDensity <- function(df, seq_col, ptm_col, int_col, format= c('wide', 'long'), save_plot= FALSE) {

  df_check <- df |> dplyr::select({{int_col}})

  if (ncol(df_check) == 0) {
    cli::cli_abort(c('x' = "The provided column names in `int_col` argument are not found in `df`.",
                     'i' = "provide correct names using `tidyselect` functions (as `starts_with`) or  in a vector `c()`."))
  }

  #normalize before?!

  df_long <- df |>
    dplyr::select({{seq_col}},{{ptm_col}},{{int_col}}) |>
    dplyr::distinct() |> #to remove duplicates
    tidyr::pivot_longer(cols = -c({{seq_col}},{{ptm_col}}),
                        names_to = "SampleName",
                        values_to = "intensity",
                        values_drop_na = TRUE
                        ) |>
    dplyr::mutate(intensity = log2(intensity),
                  sample = gsub("abundance_", "", SampleName))



  suppressMessages(print(ggplot2::ggplot(df_long,
                         ggplot2::aes(x= intensity, y= SampleName, fill = ggplot2::after_stat(x))) +
                         ggridges::geom_density_ridges_gradient() +
                         ggplot2::scale_fill_viridis_c(name = NULL,option = "C") +
                         ggplot2::labs(y= "", x = "Log2(abundance)")+
                         ggplot2::theme_minimal(base_size = 12) +
                         ggplot2::theme(legend.title = ggplot2::element_text(size=10),
                                        legend.text = ggplot2::element_text(size=10),
                                        legend.direction = "horizontal",
                                        legend.position = "top",
                                        panel.grid.major.y = ggplot2::element_blank(),
                                        panel.grid.minor.y = ggplot2::element_blank(),
                                        #axis.title.x = element_text(margin = margin(10,0,0,0), size = 12, color =  "black", face = "bold"),
                                        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(0,0,0,0), size =12, face= "bold", color = "grey60"),
                                        axis.text.x = ggplot2::element_text(vjust = 0.5, margin = ggplot2::margin(0,0,0,0), size = 12, face= "bold", color = "grey60")


                         ))) -> plot

if(save_plot){

  ggplot2::ggsave(filename =  "denisty_plot.png", plot, width= 7, height = 12, dpi = 300, bg= "white")
}


return(invisible(plot))
}
