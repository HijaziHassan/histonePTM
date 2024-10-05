





#' Calculate Coefficient of Variation
#'
#' @param df dataframe
#' @param df_meta dataframe containing at least these two columns: `SampleName` and `Condition`.
#' If `TechnicalReplicate` exist and needed to be treated separately, the easiest way to transform your `df` into long format,
#' then use `format = long` and add replace the `...` argument with the `TechnicalReplicate` column name.
#' @param seq_col sequence column
#' @param ptm_col PTM column
#' @param int_col intensity columns
#' @param format  'wide' or 'long'
#' @param ... to add any other columns
#'
#' @importFrom rlang is_missing is_empty sym
#' @importFrom dplyr select summarise filter group_by arrange left_join
#' @importFrom cli cli_abort
#' @return The input dataframe plus 3 columns per condition: `sd_condition`, `avg_condition` and `CV_condition`.
#' @export

quant_coefVariation <- function(df, df_meta, int_col, seq_col, ptm_col,  format = c("wide", "long"), ...){



  if(rlang::is_missing(df)){
    cli::cli_abort('Did you forget to pass your dataset to the `df` argument?')}
  if(rlang::is_missing(seq_col)){
    cli::cli_abort('You did not provide the name of the sequence column in `seq_col`.')}

### Wide format ---------------------
  if(format == 'wide'){  #&& save_plot== TRUE


#In wide format, metadata must be provided.
    if(rlang::is_missing(df_meta)) cli::cli_abort(c( "If the data is in `format= 'wide', the `df_meta` argument must be provided.",
                                           "x" = "`df_meta` argument is missing.",
                                           "i" = "Provide a dataframe containing `SampleNames` and `Condition` columns.
                                           SampleNames must contain the same name of the columns containing intensity values."))




    df_cv <- df


    #store intensity columns' names
    int_cols <- df |>
      dplyr::select({{int_col}}) %>%
      names()


    df_meta <- df_meta |> dplyr::filter(SampleName %in% int_cols)


    #To extract sample names from metadata as named vector
    sample_cond = split(df_meta$SampleName, df_meta$Condition)
    unique_cond = unique(df_meta$Condition)


#raise error if intensity column names are misspelled or using wrong pattern inside tidyselect function.
    if(rlang::is_empty(int_cols)) cli::cli_abort(c("Missing or incorrect input.",
                                         "x"= "The intensity columns' names were not found.",
                                         "i" = "You either forgot to add or misspelled the names you provided in the `int_col` argument." ))
    #remove NAs
    unique_cond <- unique_cond[!is.na(unique_cond)]

#iterate over each condition.
    for (cond in unique_cond) {
      #Find the indices of the intensity columns correponding to each condition
     indices = which(colnames(df_cv) %in% intersect(int_cols, sample_cond[[cond]]))

      df_cv <- df_cv |>
        mutate(

            !!paste0("avg_", cond) := rowMeans(df_cv[ ,indices, drop = FALSE], na.rm = TRUE), #drop = FALSE keeps dataframe structure otherwise it will turn to vector.
            !!paste0("sd_",  cond) := apply(df_cv[,  indices], 1, sd, na.rm = TRUE),
            !!paste0("cv_",  cond) := round((!!rlang::sym(paste0("sd_", cond)) / !!rlang::sym(paste0("avg_", cond))) * 100, 2)
        )
      }



    return(df_cv)


### Long Format -------------------------------------------
    }else if(format == "long"){



      if("Condition" %in% colnames(df)){

        df |>
          dplyr::select({{seq_col}}, {{ptm_col}}, {{int_col}}, ...) |>
          dplyr::group_by(Condition, {{ptm_col}},{{seq_col}}, ...) |>
          dplyr::summarise(CV = sd({{int_col}}, na.rm= TRUE)/mean({{int_col}}, na.rm = TRUE), .groups = "drop") |>
          dplyr::arrange(CV) -> df_cv



      }else if(rlang::is_missing(df_meta)){cli::cli_abort(c( "If the data is in `format= 'long' and contains no `Condition' column, the `df_meta` argument must be provided.",
                                                   "x" = "`Condition` column is missing.",
                                                   "i" = "Provide `df_meta` dataframe containing `SampleNames` and `Condition`."))

      }else{

        df |>
          select({{seq_col}}, {{ptm_col}}, {{int_col}}, ...) |>
          dplyr::left_join(x= _, y= df_meta, by = "SampleName") |>   #to add Condition column for coloring.
          dplyr::group_by(Condition, {{ptm_col}},{{seq_col}}, ...)
          dplyr::summarise(CV = sd({{int_col}}, na.rm= TRUE)/mean({{int_col}}, na.rm = TRUE), .groups = "drop") |>
          dplyr::arrange(CV) -> df_cv


          # tidyr::pivot_longer(cols = {{int_col}},
          #                     names_to = "SampleName",
          #                     values_to = "intensity") |>
          # dplyr::left_join(x= _, y= df_meta, by = "SampleName") |>   #to add Condition column for coloring.
          # dplyr::group_by(Condition, {{ptm_col}},{{seq_col}}, ...)
          # dplyr::summarise(CV = sd(!!!{{int_col}}, na.rm= TRUE)/mean(!!!{{int_col}}, na.rm = TRUE), .groups = "drop") |>
          # # dplyr::arrange(CV) -> df_cv


      }

      return(df_cv)

    }
#
#
#   if(plot){
#
#
#     df_long <- df_cv |>
#       pivot_longer(cols = contains(rawfilenames$SampleName),
#                    names_to = "SampleName",
#                    values_to = "intensity",
#                    values_drop_na = TRUE)
#   }
#
# #boxplot
# df_cv |>
#   ggplot(aes(x= Condition, y= CV, color = Condition))+
#
#   geom_boxplot(linewidth = .8)+
#   geom_point()+
#   geom_hline(yintercept = 0.2, linetype = "dashed", linewidth = 1.2)+
#   labs(y= "% CV", x= "")+
#   scale_color_manual(values = color_palette)+
#   scale_fill_manual(values = color_palette)+
#   scale_y_continuous(expand = c(0,0),labels = scales::percent)+
#   theme_classic(base_size = 40)+
#   coord_cartesian(clip = "off")+
#   theme(legend.position = "None",
#         axis.text.x = element_text(size= 45,colour = "black"),
#
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid.major.y = element_line(linewidth = 1, colour = "grey90"))
#
# ggsave("CV_H3_boxplot_2.png", height = 5, width = 7, dpi = 300, units =  "in", bg = "white")




  #
  #
  # #line
  #
  # df2 |>
  #   ggplot(aes(x= PTM, y= CV, color = Condition, group = Condition))+
  #   geom_point()+
  #   geom_line(linewidth = 1)+
  #   geom_hline(yintercept = 0.2, linetype = "dashed", linewidth = 1.2)+
  #   labs(y= "% CV", x= "", color= "")+
  #   facet_wrap(~variant, scales = "free")+
  #   scale_color_manual(values = color_palette)+
  #   scale_fill_manual(values = color_palette)+
  #   scale_y_continuous(labels = scales::percent)+
  #
  #   theme_classic()+
  #   theme(legend.position = "top",
  #         axis.text.y = element_text(size = 40, margin = margin(t = 0, r = 10, b = 0, l = 10)),
  #         axis.title.y = element_text(size= 42, face= "bold"),
  #         axis.text.x = element_text(size= 45, angle = 60, colour = "black",hjust = 1),
  #         legend.text = element_text(face = "bold", size= 35),
  #         strip.text = element_text(face = "bold", size= 35))
  #
  # ggsave("CV_H3.png", height = 8, width = 12, dpi = 300, units =  "in", bg = "white")



}











