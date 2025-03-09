
#' H3cano vs H3.3 Relative Abundance Plot
#' @description
#' H3K27R40 of H3cano and H3.3 are hard-coded (i.e. 'KSAPATGGVKKPHR' and 'KSAPSTGGVKKPHR').
#' These two sequences are used as proxy to the abundance of their parent protein.
#' This gives a rough estimate as some identifications are false positives.
#'
#'
#' @param df dataframe in wide format
#' @param seq_col sequence column
#' @param seq_ptm_col column name containing sequnce information
#' @param int_col column name containing peptide intensity values.
#' @param save_file \code{bool, FALSE} (default)
#' @param save_plot \code{bool, FALSE} (default)
#'
#' @importFrom dplyr select distinct filter across summarise mutate
#' @import ggplot2
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggtext element_markdown
#'
#' @return a list of two elements: the data and the plot which is a stacked bar plot showing percentage of H3.1 and H3.3 in each sample.
#' @export
#'

plot_H3vsH33 <- function(df, seq_col, seq_ptm_col, int_col, save_plot = FALSE, save_file = FALSE){

  #utils::globalVariables(c('KSAPATGGVKKPHR', 'KSAPSTGGVKKPHR'))


  RA <- df |>
    dplyr::select({{seq_col}}, {{seq_ptm_col}}, {{int_col}}) |>
    dplyr::filter({{seq_col}} %in% c('KSAPATGGVKKPHR', 'KSAPSTGGVKKPHR')) |>
    dplyr::distinct() |>

    dplyr::summarise(dplyr::across({{int_col}},
                                   \(x) sum(x, na.rm = TRUE)),
                     .by = {{seq_col}}
    ) |>
    tidyr::pivot_longer(
      cols = -{{seq_col}},
      names_to = "sample",
      values_to = "intensity",
      values_drop_na = TRUE

    )



  RA_per <- RA |>
    tidyr::pivot_wider(id_cols = sample,
                names_from = {{seq_col}},
                values_from = intensity) |>
    dplyr::mutate(per_H3_K27R40 = 100*KSAPATGGVKKPHR/(KSAPATGGVKKPHR+KSAPSTGGVKKPHR),
           per_H3.3_K27R40 = 100*KSAPSTGGVKKPHR/(KSAPATGGVKKPHR+KSAPSTGGVKKPHR),
           reorder_col = 100*KSAPATGGVKKPHR/(KSAPATGGVKKPHR+KSAPSTGGVKKPHR),
           .keep = "unused"
    )   |> tidyr::pivot_longer(cols = -c(sample, reorder_col) ,
                        names_to = "variant",
                        values_to = "intensity",
                        values_drop_na = TRUE)

  RA_per |>
    dplyr::mutate(variant= factor(variant),
                  sample = gsub('abundance_', '', sample)
    ) -> RA_per_plot

  if(save_file){

    datasave <- RA_per_plot |> dplyr::select(-reorder_col)

    openxlsx::write.xlsx(x = datasave, file = 'H3_H33.xlsx', sheetName = 'H3_H33')
  }

 plotH3 <-  RA_per_plot |>
    ggplot2::ggplot(aes(x= stats::reorder(x= sample, reorder_col, .desc = TRUE),  y= intensity, fill= variant))+
    ggplot2::geom_col(width = 0.5)+
    ggplot2::labs(x= "", y = "", fill= "",
         title = "% Relative abundance of <span style='color: #d95f02;'>H3</span> vs <span style='color: #1b9e77;'>H3.3</span>")+
    #caption= "<span style='color: #6B7787;'>Values are sum of H3K27R40 peptidoforms. Raw intensities were used.</span>")+
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::coord_flip()+
    ggplot2::geom_text(aes(label= paste(round(intensity, 2), " %")),
              color = "white", size = 5,
              position = position_stack(vjust = 0.5))+
    ggplot2::scale_y_discrete(expand = c(0,0))+
    ggplot2::theme_minimal(base_size= 12)+
    ggplot2::theme(legend.position = "none",
          axis.text.x = ggplot2::element_blank(),
          panel.grid= ggplot2::element_blank(),
          plot.title = ggtext::element_markdown(hjust = 0.5,face="bold", size = 17),
          plot.caption = ggtext::element_markdown(hjust = 0,face="italic", size = 12)
    )

if(save_plot){

ggplot2::ggsave("H3_H33.png")
}

 return(list( data= datasave , plot = plotH3))
}
