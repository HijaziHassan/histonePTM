



#' Count and plot amino aicd occurence
#' @description
#' Count the occurence of amino aicd residues in a protein sequence and plot their relative percentages.
#'
#' @param ... Protein sequence as a character string.
#' @param plot A \code{logical} argument to decide whether to plot or not the resuls.
#' @import stringr
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @import scales
#' @import tibble
#'
#' @return dataframe with 3 columns: name of a sequence, its amino acid residues and their count.
#' a plot if \code{plot} argument is set to \code{TRUE} .
#' @examples
#'H3_seq <-  "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRL
#'            VREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA"
#'H4_seq <-  "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAV
#'            TYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG"
#'seq_AAcountplot(H3_seq)
#'seq_AAcountplot(H3_seq, H4_seq, plot = TRUE)
#' @export
seq_AAcountplot <- function(..., plot = FALSE) {

seq_list <- dplyr::lst(...)

SEQsplit  <- purrr::map(unlist(seq_list), ~stringr::str_split_1(., pattern = ""))

SEQunique <- purrr::map(SEQsplit, ~stringr::str_unique(.))

df_count  <-  purrr::map2(.x = seq_list,
                          .y = SEQunique,
                          .f = ~ tibble::tibble(a.a = .y,
                                        count = round(stringr::str_count(string = .x ,
                                                                         pattern = .y)/nchar(.x), 2)*100

                          )) |> dplyr::bind_rows( .id= "seq")



if(plot == TRUE){

  plot_aa_freq <- df_count |>
    ggplot2::ggplot(ggplot2::aes(x= a.a, y= count, fill= seq)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(width  = 0.2, preserve = "single"))+
    ggplot2::labs(y= 'Relative Percentage (%)')+
    ggplot2::scale_fill_brewer(palette = "Dark2")+
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1))+
    ggplot2::theme_minimal(base_size = 14)

  print(plot_aa_freq)


}

return(df_count)

}


#'
