


#' A ridgeline visualization of intensities per sample
#' @description
#' A plot the gives an eye-bird view of the intensity distribution in each sample.
#'
#' @param df Dataset in wide format.
#' @param int_col The intensity columns. Either a vector `c('column1', 'column2', ...)` or using tidyselect functions
#' @param save_plot \code{bool}, FALSE (default).
#' like `starts_with('abundant_')` or `contains('intensity_')`
#' @importFrom ggridges geom_density_ridges_gradient
#' @import ggplot2
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom dplyr select distinct mutate everything any_of
#'
#' @return ridgeline plot
#' @export
#'


plot_intDensity <- function(df, int_col, save_plot= FALSE) {
  df_long <- df |>
    dplyr::select(dplyr::any_of({{int_col}})) |>
    dplyr::distinct() |> #to remove duplicates
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "sample", values_to = "abundance") |>
    tidyr::drop_na() |>
    dplyr::mutate(abundance = log2(abundance),
                  sample = gsub("abundance_", "", sample))

  suppressMessages(print(ggplot2::ggplot(df_long,
                         ggplot2::aes(x= abundance, y= sample, fill = ggplot2::after_stat(x))) +
                         ggridges::geom_density_ridges_gradient(size = 0.9) +
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
