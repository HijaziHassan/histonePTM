
#' Plot retention time window ranges for a PRM scheduled experiment.
#'
#' @param df A dataframe containing at least three columns: m/z and start and end retention times.
#' @param tr_start starting retention time column name.
#' @param tr_end ending retention time column name.
#' @param y_axis column name of the variable to appear on x-axis. Default is "m/z" but the axis texts are balnked out.
#' @param label Default is m/z. number will appear in the middle of the horizontal bar.
#' @param ... For any additions
#' @param save_plot \code{logical} in case needed to save the plot
#'
#' @return plot with ranges representing intervals when each precursor m/z is monitored.
#' @import ggplot2
#' @importFrom cli cli_alert_success
#'
#' @export
plot_schdueldPRMranges <- function(df,tr_start, tr_end, y_axis,label, save_plot = FALSE, ...){



p <- ggplot2::ggplot(df, ggplot2::aes(xmin = {{tr_start}},
                                xmax = {{tr_end}},
                                x= ({{tr_start}} + {{tr_end}})/2,
                                y = reorder(factor({{y_axis}}), {{tr_end}}),
                                ...

                                )
                ) +

ggplot2::geom_linerange(linewidth = 4 , alpha = 0.5, colour = "#FFD21F",
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_text(ggplot2::aes(label = if (is.numeric({{label}})) round({{label}}, 2) else {{label}}), size = 3.5, color = "black")+
    # ggplot2::scale_x_continuous(breaks = seq(0,
    #                                          floor(max(df |> dplyr::pull(tr_end))),
    #                                          5)) +
#ggplot2::scale_color_viridis_d() +
ggplot2::labs(y= "", x= " Retention time (min)")+
ggplot2::theme_classic(base_size = 14)+
ggplot2::theme(legend.title = ggplot2::element_blank(),
               #legend.position = "top",
               panel.grid.major.x = ggplot2::element_line(linetype = "dashed", colour = "grey70"),
               axis.text.y = ggplot2::element_blank(),
               axis.ticks.y = ggplot2::element_blank(),
               axis.line.y = ggplot2::element_blank(),
               axis.line.x = ggplot2::element_line(linewidth = 1.25),
                                                          ...)
 if(isTRUE(save_plot)){

    ggplot2::ggsave("prm_scheduled_ranges.png", height = 10, width =7, dpi = 300)
    cli::cli_alert_success(text = "Scheduled PRM range plot is saved successfully.")
  }

  return(p)


}



