




#' Assessing (over)Labeling Efficiency
#'
#' Sometimes overlabeling is very abundant which decrease sensitivity and affect quantification accuracy.
#'
#' @param df Proteomics results dataset.
#' @param seq_col column containing bare peptide sequences
#' @param seq sequence(s) to base the assessment on. Multiple sequences should be passed as a character vector.
#' If empty or `NULL` (default), five H3 and H4 built-in sequences will be analyzed.
#' If "All", then all the unique sequence in the `seq_col` will be analyzed.
#' @param ptm_col column containing peptide ptm information (for overTMA: TMAyl, tmaNt, and tmame1, tma;, for overProp: Propionyl, prNt, pr). Only
#'two representations of PTMs can be recognized to filter upon: Propionyl (K12) or K27bu-S32pr ....)
#' @param int_col sample columns containing intensity values
#' @param save_plot bool; `FALSE` (default). Save or not the plot to the desktop.`
#' @param plot_title A title for the plot
#' @importFrom stringr str_detect str_count str_remove
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate summarise across select distinct if_else cur_column pull bind_rows case_when
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom cli cli_alert_info cli_alert_warning
#' @import ggplot2
#'
#' @return list of two dataframes (raw and normalized values) and a bar plot.
#' @export


ptm_labelingAssessment <- function(df, seq_col, seq = NULL, ptm_col, int_col, save_plot = FALSE, plot_title= ""){



  #STY propionylation/TMAylation
  regex_OverLab = 'Propionyl.*\\([TSY]\\d+\\)|TMAyl(?:_correct)?.*\\([TSY]\\d+\\)|[ST]\\d+(?:tma|pr)'

  #Methyl or me1 and not tmame1
  regex_nonLabMe1 = '(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:[:upper:]*\\d+me1)\\b'
    #^(?:(?!.*(?:\\(Any N-term\\)|prNt|tmaNt)).)*$|(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:[:upper:]*\\d+me1)\\b

  #how many Ks are modified
  regex_Kmod = '\\(K\\d+?\\)|K\\d+'

  #Propionyl//TMA N-term
  regex_Nterm = '.* \\(Any N-term\\)|tmaNt|prNt'



    if (identical(seq, "All")) {
      seq <- unique(df[[as.character(substitute(seq_col))]])
    } else if (is.null(seq) || length(seq) == 0) {
      seq <- c('ISGLIYEETR', 'YRPGTVALREIR', 'KQLATKAAR', 'DNIQGITKPAIR', 'DAVTYTEHAKR', 'KSTGGKAPR', 'KSAP.TGGVKKPHR')
    }

    df_filtered <- df |>
      dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}) |>
      dplyr::filter(stringr::str_detect({{seq_col}}, paste(seq, collapse = "|"))) |>
      dplyr::distinct(dplyr::across({{int_col}}), .keep_all = TRUE)


    if (nrow(df_filtered) == 0) {
      cli::cli_alert_warning("None of the specified sequences were found in the data.")
      return(NULL)
    }



    df_results <- purrr::map(seq, function(current_seq) {
      df_seq <- df_filtered |>
        dplyr::filter(stringr::str_detect({{seq_col}}, current_seq))

      if (nrow(df_seq) < 1) {
        cli::cli_alert_info('The sequence {.val {current_seq}} is not found in your {.arg seq_col} column.')
        return(tibble::tibble())
      }

      df_tagged <- df_seq |>
        dplyr::mutate(
          isOverLab = stringr::str_detect({{ptm_col}}, regex_OverLab),
          isNonLabeledme1 = stringr::str_detect({{ptm_col}}, regex_nonLabMe1),
          isFullyModified = dplyr::if_else(
            stringr::str_count({{seq_col}}, "K") == stringr::str_count({{ptm_col}}, regex_Kmod),
            TRUE, FALSE
          ),
          isNterm = stringr::str_detect({{ptm_col}}, regex_Nterm)
        )


      df_total <- df_tagged |>
        dplyr::summarise(dplyr::across({{int_col}}, ~sum(.x, na.rm = TRUE)))



      df_all <- df_tagged |>
        dplyr::mutate(
          label = dplyr::case_when(
            isNonLabeledme1 | !isFullyModified | !isNterm & !isOverLab ~ 'UnderLabeled',
            isOverLab & !(isNonLabeledme1 | !isFullyModified | !isNterm) ~ 'OverLabeled',
            isNterm & !isOverLab & !isNonLabeledme1 & isFullyModified ~ 'Desired',
            .default = 'UnderOverLabeled'
          )
        ) |>

        dplyr::summarise(dplyr::across({{int_col}}, ~sum(.x, na.rm = TRUE)), .by = label) |>
        dplyr::mutate(seq_analyzed = current_seq)

      df_all_norm <- df_all |>
        dplyr::mutate(
          dplyr::across(
            .cols = {{int_col}},
            .fns = ~ . * 100 / df_total |> dplyr::pull(dplyr::cur_column())
          )
        )



      # Plot for each sequence

      plot <- df_all_norm |>
        dplyr::mutate(label = factor(label, levels= c('Total', 'Desired', 'OverLabeled', 'UnderLabeled', 'UnderOverLabeled'))) |>
        tidyr::pivot_longer(
          cols = {{int_col}},
          names_to = 'sample',
          values_to = 'intensity'
        ) |>
        dplyr::mutate(sample = stringr::str_remove(sample, 'abundance_')) |> # specific for Proline software
        ggplot2::ggplot(ggplot2::aes(x= sample, y= intensity, fill= label)) +
        ggplot2::labs(y= "Relative Intensity (%)", fill= "", x= "", title= plot_title) +
        ggplot2::geom_col(position= ggplot2::position_dodge2(preserve= 'single')) +
        ggplot2::geom_hline(yintercept = 100,  linetype = "dotted", color= "#DEDEDE")+
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text = ggplot2::element_text(face= "bold"),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5),
          axis.line = ggplot2::element_line(linewidth = 1.2),
          axis.text = ggplot2::element_text(face= 'bold', color = 'black')
        )

      if(save_plot){
      cat(paste0('Plotting: ', current_seq))
      filename <- paste0('overlabassess_', current_seq, ".png")
      ggplot2::ggsave(filename, plot, width = 8, height = 6)
      }
      return(df_all_norm)

    }) |> dplyr::bind_rows()

    return(df_results)
  }



