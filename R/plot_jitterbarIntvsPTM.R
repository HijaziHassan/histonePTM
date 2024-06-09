#' Intensity vs PTM grouped bar plot with jittered points
#' @description Each individual measurement will be plotted as datapoint.
#' Their median/mean will be represented as a bar.
#' Each bar will be colored according to the condition.
#'
#' @param dataset your data loaded as dataframe. Must be in long format.
#' @param x_axis x variable (PTM column)
#' @param y_axis y variable (intensity column). If already in %, set \code{scale} to 1.
#' @param condition The condition column (WT vs disease, concentration, ...)
#' @param sequence peptide sequence column
#' @param label peptide label to be used for y-axis title (e.g. 'H3.3')
#' @param title_plot A plot title (optional)
#' @param fun median (default) or mean. This will be the height of the bar.
#' @param save_plot (\code{logical})
#' @param scale 100 (default). If you want to keep values as they are use 1.
#'
#' @importFrom stats reorder
#' @importFrom scales label_percent
#' @import ggplot2
#'
#' @return ggplot object
#' @export
#'
plot_jitterbarIntvsPTM <- function(dataset,
                             x_axis, #PTM
                             y_axis, #Intensity
                             condition, #variable on which comparison is done
                             sequence, # peptide sequence col
                             label= "", # labels will be used for y-axis title
                             plot_title ="", #optional
                             fun = c("median", "mean"),
                             scale = 100,
                             save_plot = FALSE
){

  fun = match.arg(fun)



  p <- ggplot2::ggplot(dataset,
                  ggplot2::aes(    x = stats::reorder({{x_axis}},{{y_axis}}),
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
                         shape= 21, color = "black", size = 1.1, alpha = 1 ) +

    ggplot2::scale_y_continuous(labels = scales::label_percent(scale = scale)) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set1") +

    ggplot2::labs(x= "",
                  y= paste0('% of all variably modified forms of ', label),
                  title = plot_title) +


    ggplot2::theme(
         panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.grid.major.y  = ggplot2::element_line(color = "grey70"),
              axis.text.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 20)),
             axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
              axis.text.x = ggplot2::element_text(size = 12, angle = 60, colour = "black", hjust = 0.7),
              axis.line.x = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              plot.margin = ggplot2::unit(c(0.5, 0.5, 0.2, 0.5), "cm") ,
              legend.text = ggplot2::element_text(face = "bold", size = 10),
              legend.position.inside = c(0.01, 0.8),
                legend.justification = c(0, 0),
                    legend.direction = "vertical",
            legend.title = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
                   title = ggplot2::element_text(size= 16)
    )


   writeLines(paste0("Plotting: ", label))


   if(save_plot){
     ggplot2::ggsave(filename = paste0(label, ".png"),
            dpi = 300,
            height = 7,
            width = 10,
            units =  "in",
            bg = "white")
   }

return(p)
}
