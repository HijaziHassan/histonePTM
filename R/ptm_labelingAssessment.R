




#' Assessing (over)Labeling Efficiency
#'
#' Derivatization efficiency needs to be checked. Sometimes over/underlabeling is very abundant which decreases sensitivity and affects quantification accuracy.
#'
#' @param df Proteomics results dataset. STYlabel must be set as variable modification in addition to K(me1+label).
#' @param seq_col column containing bare peptide sequences
#' @param seq sequence(s) to base the assessment on. Multiple sequences should be passed as a character vector.
#' If empty or  \code{NULL} (default), five H3 and H4 built-in sequences will be analyzed.
#' If "All", then all the unique sequence in the \code{seq_col} will be analyzed.
#' @param ptm_col column containing peptide ptm information (for overTMA: TMAyl, tmaNt, and tmame1, tma; for overProp: Propionyl, prNt, pr). Only
#'two representations of PTMs can be recognized to filter upon: Propionyl (K12) or K27bu-S32pr ....)
#' @param int_col sample columns containing intensity values.
#' @param type Either "dodged (default)" or "stacked" to plot stacked or dodged bar plot.
#' @param show_text bool;  \code{TRUE} (default) shows the percentages on the bars.
#' The values are rounded to the nearest whole number unless it less than one, then it is rounded to the nearest tenth.
#' @param save_plot bool;  \code{FALSE} (default). Save or not the plot(s).`
#' @param plot_title A title for the plot
#' @param output_dir character; A folder name to store the plots into if \code{save_plot} is set \code{TRUE}. If \code{NULL} (default), plots will be saved into the working directory.
#' @importFrom stringr str_detect str_count str_remove
#' @importFrom tidyr pivot_longer complete
#' @importFrom dplyr filter mutate summarise across select distinct if_else cur_column pull bind_rows case_when all_of
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom cli cli_alert_info cli_alert_warning cli_inform
#' @import ggplot2
#'
#' @return A list of three: a dataframe combining all the data and a list of plot(s) for each sequence, and a dataframe for ID counts.
#' @details
#' Labeling categories are defined as follows:
#'
#' \describe{
#'   \item{\strong{Desired}}{Peptides in which all lysine residues are modified or labeled in addition to the N-term.}
#'   \item{\strong{OverLabeled}}{Peptides that have at least one labeled serine (S), threonine (T), or tyrosine (Y) residue.}
#'   \item{\strong{UnderLabeled}}{
#'     Peptides that contain any of the following unlabeled features:
#'     \itemize{
#'       \item Unlabeled peptide N-terminus
#'       \item Unlabeled monomethylated lysine (Kme1)
#'       \item Unlabeled lysine (K)
#'     }
#'     These may occur individually or together.
#'   }
#'   \item{\strong{UnderOverLabeled}}{Peptides that meet the criteria for both `UnderLabeled` and `OverLabeled`.}
#' }


#'
#' @export


ptm_labelingAssessment <- function(df,
                                   seq_col,
                                   seq = NULL,
                                   ptm_col,
                                   int_col,
                                   type = c('dodged', 'stacked'),
                                   show_text = TRUE,
                                   save_plot = FALSE,
                                   plot_title= "",
                                   output_dir= NULL){

type = match.arg(type)


#store intensity columns' names
int_cols <- df |>
  dplyr::select(dplyr::all_of({{int_col}}))  |>
  names()

#raise error if intensity column names are misspelled or using wrong pattern inside tidyselect function.
if(rlang::is_empty(int_cols)) cli::cli_abort(c("Missing or incorrect input.",
                                               "x"= "The intensity columns' names were not found.",
                                               "i" = "You either forgot to add or misspelled the names you provided in the `int_col` argument." ))

  #STY propionylation/TMAylation
  regex_OverLab = 'Propionyl.*\\([TSY]\\d+\\)|TMAyl(?:_correct)?.*\\([TSY]\\d+\\)|[TSY]\\d+(?:tma|pr)'

  #Methyl or me1 and not tmame1
  regex_nonLabMe1 = '(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:K\\d+me1)\\b'
    #^(?:(?!.*(?:\\(Any N-term\\)|prNt|tmaNt)).)*$|(?:Methyl\\s*K?R?\\s*\\(K\\d+\\))|(?:[:upper:]*\\d+me1)\\b

  #how many Ks are modified
  regex_Kmod = '\\(K\\d+?\\)|K\\d+'

  #Propionyl//TMA N-term
  regex_Nterm = '.* \\(Any N-term\\)|.* \\(Protein N-term\\)|tmaNt|prNt|acPNt'



    if (identical(seq, "All")) {
      seq <- unique(df[[as.character(substitute(seq_col))]])
    } else if (is.null(seq) || length(seq) == 0) {
      seq <- c('ISGLIYEETR', 'YRPGTVALREIR', 'KQLATKAAR', 'DNIQGITKPAIR', 'DAVTYTEHAKR', 'KSTGGKAPR', 'KSAP.TGGVKKPHR')
    }

  #df containing sequence (containing all selected seqs), ptm and intensity columns
    df_filtered <- df |>
      dplyr::select({{seq_col}}, {{ptm_col}}, dplyr::all_of({{int_col}})) |>
      dplyr::filter(stringr::str_detect({{seq_col}}, paste0("^(", paste(seq, collapse = "|"), ")$"))) |>
      dplyr::distinct(dplyr::across(dplyr::all_of({{int_col}})), .keep_all = TRUE)


    if (nrow(df_filtered) == 0) {
      if(length(seq)>1){
      cli::cli_alert_warning("None of the specified sequences {.val {seq}} were found in your data.")
      }else{

        cli::cli_alert_warning("The sequence {.val {seq}} was not found in your data.")

      }

      return()
    }

# Analysis ----------------------------------

plot_list = list()
df_list <- list()
#Analyze for each sequence
    results <- purrr::map(seq, function(current_seq) {
     #dataframe containing only one sequence
       df_seq <- df_filtered |>
        dplyr::filter(stringr::str_detect({{seq_col}}, current_seq))
       #if the sequence is not found, return empty tibble
      if (nrow(df_seq) < 1) {
        cli::cli_alert_info('The sequence {.val {current_seq}} is not found in your {.arg seq_col} column.')
        return(tibble::tibble())
      }
      #analyze the PTMs adding tags columns
      df_tagged <- df_seq |>
        dplyr::mutate(
          isOverLab = stringr::str_detect({{ptm_col}}, regex_OverLab),
          isNonLabeledme1 = stringr::str_detect({{ptm_col}}, regex_nonLabMe1),
          isFullyModified = dplyr::if_else(
            stringr::str_count({{seq_col}}, "K") == stringr::str_count({{ptm_col}}, regex_Kmod),
            TRUE, FALSE
          ),
          isNterm = stringr::str_detect({{ptm_col}}, regex_Nterm),
          isFullyUnmod = is.na({{ptm_col}})
        )
# Intensity-based analysis -------------------------------

# dataframe for total intensity of each sample
      seq_total_sum <- df_tagged |>
        dplyr::summarise(dplyr::across(dplyr::all_of({{int_col}}), ~sum(.x, na.rm = TRUE)))


#based on the tag,define the labeling status
      df_labeled <- df_tagged |>
        dplyr::mutate(
          tag = dplyr::case_when(
            isNonLabeledme1 | !isFullyModified | !isNterm & !isOverLab | isFullyUnmod ~ 'UnderLabeled',
            isOverLab & !(isNonLabeledme1 | !isFullyModified | !isNterm) ~ 'OverLabeled',
            isNterm & !isOverLab & !isNonLabeledme1 & isFullyModified ~ 'Desired',
            .default = 'UnderOverLabeled'
          )
        )

      #dataframe total intensity of each category in each sample
      df_summary <- df_labeled |>
        dplyr::summarise(dplyr::across(dplyr::all_of({{int_col}}), ~sum(.x, na.rm = TRUE)), .by = tag) |>
        dplyr::mutate(seq_analyzed = current_seq)  # Add sequence column

      #normalize by total intensity of each sample
      df_all_norm <- df_summary|>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::all_of({{int_col}}),
            .fns = ~ . * 100 / seq_total_sum |> dplyr::pull(dplyr::cur_column())
          )
        )



# Counting-based analysis -------------------------------
      df_count <- df_labeled |>
        dplyr::mutate(tag = factor(tag, levels= c('Desired', 'OverLabeled', 'UnderLabeled', 'UnderOverLabeled'))) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of({{int_col}}),
          names_to = 'sample',
          values_to = 'intensity',
          values_drop_na = TRUE
        ) |>
        #dplyr::mutate(seq_analyzed = current_seq) |>
        dplyr::count({{seq_col}}, sample, tag, name = "tag_count") |>
        tidyr::complete({{seq_col}}, sample, tag, fill = list(tag_count = 0))


# Plotting -------------------------------
#dataframe with columns: tag, seq_analyzed, sample, intensity, and pct_formatted.
      plot <- df_all_norm |>
        dplyr::mutate(tag = factor(tag, levels= c('Desired', 'OverLabeled', 'UnderLabeled', 'UnderOverLabeled'))) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of({{int_col}}),
          names_to = 'sample',
          values_to = 'intensity'
        ) |>   tidyr::complete(
          seq_analyzed,    # <— include it here
          sample,
          tag,
          fill = list(intensity = 0)
        ) |>
        dplyr::mutate(sample = stringr::str_remove(sample, 'abundance_'), # specific for Proline software
                      pct_formatted = {
                        if (!show_text) {
                          ''
                        } else if (type == "dodged") {
                          ifelse(intensity== 0, "",
                            ifelse(intensity < 1, "*",
                                 #paste0(round(intensity, digits = 1), "%"),
                                 paste0(round(intensity), "%"))
                                 )
                        } else {
                          ifelse(intensity < 5 | is.na(intensity), "", paste0(round(intensity), "%"))
                        }
                      }

                      )



      sample_number <- length(colnames(df_filtered))-2

      angle = ifelse(sample_number < 10, 0, 45)

      custom_theme <-  ggplot2::theme_classic(base_size = 11)+
        ggplot2::theme(
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text = ggplot2::element_text(face= "bold"),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.key.width = unit(0.2, 'cm'),
          legend.justification = c(0, 0),
          axis.text.x = ggplot2::element_text(angle = angle, vjust = 0.5),
          axis.title.y = ggplot2::element_text(face = "bold", size = 12),
          axis.line = ggplot2::element_line(linewidth = 1.2),
          axis.text = ggplot2::element_text(face= 'bold', color = 'black')
        )


      barplotting <- function(){

        ggplot2::ggplot(plot, ggplot2::aes(x= sample, y= intensity, fill= tag, label= pct_formatted)) +
          ggplot2::labs(y= "Relative Intensity (%)", fill= "", x= "", title= plot_title) +
          ggplot2::geom_col(position= ggplot2::position_dodge2(preserve= 'single', width = 0.8)) +
          ggplot2::geom_text(size= 3.5, fontface = 'bold', color= 'black',
                             , position = ggplot2::position_dodge(0.9), vjust= -.5)+
          ggplot2::geom_hline(yintercept = 100,  linetype = "dotted", color= "#DEDEDE")+
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +
          #ggplot2::scale_fill_brewer(palette = "Dark2") +
          ggplot2::scale_fill_manual(values= c("Desired" = "#1B9E77",
                                       "OverLabeled" = "#D95F02",
                                       "UnderLabeled" = "#7570B3",
                                       "UnderOverLabeled" = "#E7298A"))+
          ggplot2::coord_cartesian(clip = 'off')+
          custom_theme

      }


      stackplotting <- function(){

        ggplot2::ggplot(plot, ggplot2::aes(x= sample, y= intensity, fill= tag, label= pct_formatted)) +
          ggplot2::labs(y= "Relative Intensity (%)", fill= "", x= "", title= plot_title) +
          ggplot2::geom_col(width = 0.6) +
          ggplot2::geom_text(position= ggplot2::position_stack(vjust= 0.5),
                             size= 4, fontface = 'bold', color= 'white', na.rm = TRUE)+
          ggplot2::geom_hline(yintercept = 100,  linetype = "dotted", color= '#b4b4b4')+
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +
          ggplot2::scale_fill_manual(values = c("Desired" = "#1B9E77",
                                       "OverLabeled" = "#D95F02",
                                       "UnderLabeled" = "#7570B3",
                                       "UnderOverLabeled" = "#E7298A"))+
          custom_theme
      }


      plot_functions <- list(
        dodged = barplotting,
        stacked = stackplotting
      )


      plot_output <- switch(
        type,
        dodged = barplotting(),
        stacked = stackplotting(),
        stop("Invalid plot type. Use 'dodged' or 'stacked'.")
      )



      if(save_plot){
        output_dir <- if (is.null(output_dir)) getwd() else output_dir
        if (!is.null(output_dir) && !dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
          cli::cli_inform('The folder {.val {output_dir}} is created to store the plot(s).')
        }
        cli::cli_inform("Plotting: {.val {current_seq}}")
        filename <- paste0('overlabassess_', current_seq, ".png")
        ggplot2::ggsave(filename = filename, plot = plot_output,path = output_dir, width = 8, height = 6, dpi = 300)
      }



      list(
        df_all_norm = df_all_norm,
        plot_output = plot_output,
        df_count= df_count
      )

    }


    )

    names(results) <- seq

    # Extract all dataframes and bind them together
    df_all_norm_combined <- results |>
      purrr::map("df_all_norm") |>
      dplyr::bind_rows()

    # Extract all plots into a list
    plot_list <- results |>
      purrr::map("plot_output")

    df_count <- results |>
      purrr::map("df_count")


    return(list(
      data = df_all_norm_combined,
      plot = plot_list,
      counting = df_count
    ))


}




