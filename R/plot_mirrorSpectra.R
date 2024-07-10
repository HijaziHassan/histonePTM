




.read_any_format <- function(file) {
  # Check file extension
  file_ext <- tolower(tools::file_ext(file))

  if (file_ext %in% c("csv", "txt")) {

    data <- read.table(file, header = TRUE, sep = ",")
  } else if (file_ext %in% c("xls", "xlsx")) {

    data <- readxl::read_excel(file)
  } else {

    cli::cli_abort(c("Unsupported file format.",
                     "i"= "The {file_ext} cannot be read properly." ,
                     "x"= "Please provide a csv, txt, or excel file."))
  }

  return(data)
}



#' Mirror Plot
#'
#' @param data A \code{dataframe} or a comma-separated (txt, csv) or excel (.xls, xlsx) file.
#' @param mz_col mz values column
#' @param intensity_col intensity values column
#' @param grouping_col Any column that serve as ID. Could be scan number, or file, ... Those will annotate the y-axis.
#' @param top_spec Value (\code{character}) from the \code{grouping col}. Ex: 'processed_spec".
#' @param bottom_spec Value (\code{character}) from the \code{grouping col}. Ex: 'raw_spec".
#' @param title Title of the plot.
#' @param norm_basepeak \code{logical}. If \code{TRUE} (default), intensities will be normalized to base peak intensity.
#' @param readfromfile \code{logical}. If \code{FALSE} (default), it means that the \code{data} argument is  a \code{dataframe} and
#' not a file name (\code{.csv}, \code{txt}, \code{xls} or \code{xlsx}.
#'
#' @import stringr
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import readxl
#' @import cli
#'
#' @return A mirror plot.
#' @export
plot_mirrorSpectra <- function(data,
                               mz_col, intensity_col, grouping_col,
                               top_spec , bottom_spec, title= "",
                               norm_basepeak =TRUE,  readfromfile = FALSE){

  #read from file
  if(isTRUE(readfromfile)){
    spectra <- .read_any_format(data) |> tibble::as_tibble()
  }else{
    spectra <- data
  }

  #read from R object
  spectra <- spectra |>
    dplyr::mutate({{intensity_col}} := as.numeric({{intensity_col}}),
           {{mz_col}} := as.numeric({{mz_col}}),
           {{grouping_col}} := as.character({{grouping_col}})
           ) |>
    dplyr::filter(!{{intensity_col}} == 0)   #remove m/z (i.e. peaks) with 0 intensity

  #normalize to base peak
  if(isTRUE(norm_basepeak)){
    spectra <- spectra |>
      dplyr::mutate({{intensity_col}} := {{intensity_col}} / max({{intensity_col}}) * 100, .by = {{grouping_col}})
  }


  #reverse bottom spectrum
  spectra <- spectra |>
    dplyr::mutate({{intensity_col}} := ifelse({{grouping_col}} == {{bottom_spec}} ,-intensity, intensity))


#plot
  pal <- c("#d1495b", "#00b159")
  selected_spec <- c({{bottom_spec}}, {{top_spec}})
  names(pal) <-   selected_spec


  y_title = stringr::str_glue("<span style = 'color:#ffffff;'>.........</span>","<b style='color:#d1495b'>{bottom_spec}</b>",
                     "<span style = 'color:#ffffff;'>.................</span>",
                     "<b style='color:#00b159'>{top_spec}</b>")


  mirroplot <- ggplot2::ggplot(spectra  |>
           dplyr::filter( {{grouping_col}} %in%  selected_spec), ggplot2::aes({{mz_col}}, {{intensity_col}},
                                                              xend = {{mz_col}},
                                                              yend = 0,
                                                              color = {{grouping_col}}))+
    ggplot2::geom_segment(linewidth = 1) +
    ggplot2::labs(y = y_title, title = {{title}})+
    ggplot2::scale_color_manual(values = pal)+
    ggplot2::theme_bw(base_size =  12) +
    ggplot2::theme(
      legend.position = "none",
           plot.title = ggplot2::element_text(hjust = 0.5),
         axis.title.x = ggplot2::element_blank(),
         axis.ticks.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
         axis.title.y = ggtext::element_markdown()
    )



return(mirroplot)
}









