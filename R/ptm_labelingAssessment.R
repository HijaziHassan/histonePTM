




#' Assessing (over)Labeling Effeciency
#'
#' Some time overlabeling is very abundant which decrease sensitivity and affect quantification accuracy.
#'
#' @param df Your dataset
#' @param seq_col column containing bare peptide sequences
#' @param seq sequence to base the assessment on it
#' @param ptm_col column containing peptide ptm information
#' @param int_col sample columns containing intensity values
#' @param plot_title A title for the plot
#' @importFrom stringr str_detect str_count str_remove
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @import ggplot2
#'
#' @return list of two dataframes (raw and normalized values) and a bar plot.
#' @export


ptm_labelingAssessment <- function(df, seq_col, seq, ptm_col, int_col, plot_title= ""){



  #STY propionylation/TMAylation
  regex_OverProp = 'Propionyl.*\\([TSY]\\d+\\)|TMAyl_correct.*\\([TSY]\\d+\\)'

  #No N-term prop or methyl
  regex_UnderProp = '^(?:(?!.*\\(Any N-term\\)).)*$|(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))'

  #how many Ks are modified
  regex_Kmod = '\\(K\\d+\\)'

  #Propionyl//TMA N-term
  regex_Nterm = '.* \\(Any N-term\\)'

  df_filtered <- df |>
    dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}) |>
    dplyr::filter(stringr::str_detect({{seq_col}}, {{seq}})) |>
    #filter(complete.cases(across(starts_with('abundance_')))) |>
    dplyr::distinct(dplyr::across({{int_col}}), .keep_all = TRUE)







  df_tagged <- df_filtered |>
    dplyr::mutate(isOverProp = stringr::str_detect({{ptm_col}}, regex_OverProp),
                  isUnderProp = stringr::str_detect({{ptm_col}}, regex_UnderProp),
                  isFullyModified= dplyr::if_else(stringr::str_count({{seq_col}}, "K") == stringr::str_count({{ptm_col}}, regex_Kmod),
                                                  "TRUE",
                                                  "FALSE"),
                  isNterm = stringr::str_detect({{ptm_col}}, regex_Nterm ))


  df_total <- df_tagged |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'Total')



  df_underprop <- df_tagged |>
    dplyr::filter(isUnderProp == TRUE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'UnderLabeled')

  df_overprop <- df_tagged |>
    dplyr::filter(isOverProp == TRUE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'OverLabeled')


  df_desired <- df_tagged |>
    dplyr::filter(isNterm == TRUE , isFullyModified == TRUE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'Desired')



  df_total = dplyr::bind_rows(df_total, df_underprop, df_overprop, df_desired)


  df_total_norm <- df_total |>
    dplyr::mutate( dplyr::across(.cols =  {{int_col}},
                                 .fns = ~ .*100 / df_total |>
                                   dplyr::filter(label == 'Total') |>
                                   dplyr::pull(dplyr::cur_column())
    )
    )


  #Change angle of x-axis labels according to sample number

  sample_number <- length(colnames(df_filtered))-2

  angle = ifelse(sample_number < 4, 0, 45)


  p <- df_total_norm |>
    dplyr::mutate(label = factor(label, levels= c('Total', 'Desired', 'OverLabeled', 'UnderLabeled'))) |>
    tidyr::pivot_longer( cols = {{int_col}},
                         #names_pattern = 'abundance_(.*)',
                         names_to = 'sample',
                         values_to = 'intensity' ) |>
    dplyr::mutate(sample= stringr::str_remove(sample, 'abundance_')) |>
    ggplot2::ggplot(aes(x= sample, y= intensity, fill= label))+
    ggplot2::labs(y= "% Relative Intensity", fill= "", x= "", title= plot_title)+
    ggplot2::geom_col(position= ggplot2::position_dodge2(preserve= 'single'))+
    ggplot2::scale_fill_brewer(palette = "Dark2")+
    ggplot2::theme_classic(base_size = 14)+
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(face= "bold"),
                   axis.text.x = ggplot2::element_text(angle = angle, vjust = 0.5))


  #what about misclevages
  #make sure how underprop and overprop are related.
  #adapt it to take more than one sequence.

  return(list(raw_data = df_total, normalized_data = df_total_norm, plot= p))
}

