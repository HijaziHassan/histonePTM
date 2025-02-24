




#' Assessing (over)Labeling Effeciency
#'
#' Some time overlabeling is very abundant which decrease sensitivity and affect quantification accuracy.
#'
#' @param df Your dataset
#' @param seq_col column containing bare peptide sequences
#' @param seq sequence to base the assessment on it
#' @param ptm_col column containing peptide ptm information (for overTMA: TMAyl, tmaNt, and tmame1, tma;, for overProp: Propionyl, prNt, pr). Only
#'two representations of PTMs can be recognized to filter upon: Propionly (K12) or K27bu-S32pr ....)
#' @param int_col sample columns containing intensity values
#' @param plot_title A title for the plot
#' @importFrom stringr str_detect str_count str_remove
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate summarise across select distinct if_else cur_column pull bind_rows
#' @import ggplot2
#'
#' @return list of two dataframes (raw and normalized values) and a bar plot.
#' @export


ptm_labelingAssessment <- function(df, seq_col, seq, ptm_col, int_col, plot_title= ""){



  #STY propionylation/TMAylation
  regex_OverLab = 'Propionyl.*\\([TSY]\\d+\\)|TMAyl(?:_correct)?.*\\([TSY]\\d+\\)|[ST]\\d+(?:tma|pr)'

  #Methyl or me1 and not tmame1
  regex_nonLabMe1 = '(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:[:upper:]*\\d+me1)\\b'
    #^(?:(?!.*(?:\\(Any N-term\\)|prNt|tmaNt)).)*$|(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:[:upper:]*\\d+me1)\\b

  #how many Ks are modified
  regex_Kmod = '\\(K\\d+?\\)|K\\d+'

  #Propionyl//TMA N-term
  regex_Nterm = '.* \\(Any N-term\\)|tmaNt|prNt'



  # filter data by removing duplications, selecting the sequence of interest
  df_filtered <- df |>
    dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}) |>
    dplyr::filter(stringr::str_detect({{seq_col}}, {{seq}})) |>
    dplyr::distinct(dplyr::across({{int_col}}), .keep_all = TRUE)

#### Calculate CVs #####

if(nrow(df_filtered)>1){


#### Label the data #####

  df_tagged <- df_filtered |>
    dplyr::mutate(isOverLab = stringr::str_detect({{ptm_col}}, regex_OverLab),
                  isNonLabeledme1 = stringr::str_detect({{ptm_col}}, regex_nonLabMe1),
                  isFullyModified= dplyr::if_else(stringr::str_count({{seq_col}}, "K") == stringr::str_count({{ptm_col}}, regex_Kmod),
                                                  "TRUE",
                                                  "FALSE"),
                  isNterm = stringr::str_detect({{ptm_col}}, regex_Nterm ))


##### Isolate and sum each labeled group #####
  df_total <- df_tagged |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'Total')



  df_underlab <- df_tagged |>
                #monomethyl unlabeled OR unlabeled K OR not N-term modified
    dplyr::filter(isNonLabeledme1 == TRUE | isFullyModified == FALSE | isNterm == FALSE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'UnderLabeled')

  df_overlab <- df_tagged |>
    dplyr::filter(isOverLab == TRUE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'OverLabeled')


  df_desired <- df_tagged |>
    dplyr::filter(isNterm == TRUE ,
                  isOverLab == FALSE,
                  isNonLabeledme1 == FALSE,
                  isFullyModified == TRUE) |>
    dplyr::summarise(dplyr::across({{int_col}},  ~sum(.x, na.rm = T))) |>
    dplyr::mutate(label = 'Desired')



  df_total = dplyr::bind_rows(df_total, df_underlab, df_overlab, df_desired) |>
    mutate(seq_analyzed = seq)


  df_total_norm <- df_total |>
    dplyr::mutate( dplyr::across(.cols =  {{int_col}},
                                 .fns = ~ .*100 / df_total |>
                                   dplyr::filter(label == 'Total') |>
                                   dplyr::pull(dplyr::cur_column())
    )
    )


  #Change angle of x-axis labels according to sample number

  sample_number <- length(colnames(df_filtered))-2

  angle = ifelse(sample_number < 8, 0, 45)


  p <- df_total_norm |>
    dplyr::mutate(label = factor(label, levels= c('Total', 'Desired', 'OverLabeled', 'UnderLabeled'))) |>
    tidyr::pivot_longer( cols = {{int_col}},
                         #names_pattern = 'abundance_(.*)',
                         names_to = 'sample',
                         values_to = 'intensity' ) |>
    dplyr::mutate(sample= stringr::str_remove(sample, 'abundance_')) |>
    ggplot2::ggplot(aes(x= sample, y= intensity, fill= label))+
    ggplot2::labs(y= "Relative Intensity (%)", fill= "", x= "", title= plot_title)+
    ggplot2::geom_col(position= ggplot2::position_dodge2(preserve= 'single'))+
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01)))+
    ggplot2::scale_fill_brewer(palette = "Dark2")+
    ggplot2::theme_classic(base_size = 14)+
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(face= "bold"),
                   axis.text.x = ggplot2::element_text(angle = angle, vjust = 0.5),
                   axis.line = ggplot2::element_line(linewidth = 1.2),
                   axis.text = element_text(face= 'bold', color = 'black')
                   )


  #what about misclevages
  #make sure how underprop and overprop are related.
  #adapt it to take more than one sequence.

  return(list(raw_data = df_total, normalized_data = df_total_norm, plot= p))

    }else if (nrow(df_filtered) == 0){

  cli::cli_alert_info('The sequence {.val {seq}} is not found in your {.val {deparse(substitute(seq_col))}} column.')

      }else{

         return(NULL)

  }

}

