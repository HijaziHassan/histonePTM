#' A quick literature review on NCBI
#'
#' @description
#' A wrapper function around \code{rentrez} package functions to extract article titles using year and search term.
#' Check \code{rentrez package} \href{https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html}{documentation}.
#'
#' @param start An \code{integer} specifying the year from which to start the search.
#' @param end An \code{integer} specifying the year at which to stop the search.
#' @param term any search term including boolean operators (i.e. AND, NOT,OR, ...) and tags ([au] for author,
#' [lau] for last author, review [pt] for publication type, [ti] for title etc...Check \href{https://pubmed.ncbi.nlm.nih.gov/help/#using-search-field-tags}{here}
#' for an exhaustive list.
#' @param db database to searhc in. defaults is "pubmed". check \code{entrez_dbs()} for all options.
#' @param save_plot \code{bool}. If `TRUE` a line plot of count vs year colored by search term will be exported as `png`.
#' @param save_file \code{bool}. If `TRUE` results will be saved as `.csv` files.
#' @importFrom stringr str_glue str_remove_all
#' @importFrom purrr map map_chr map_dbl
#' @importFrom dplyr mutate select summarize filter group_by across where
#' @importFrom tidyr tibble unnest
#' @importFrom utils write.csv
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom scales pretty_breaks
#' @import ggplot2
#'
#' @return A list. A dataframe with  year, id, and article title. A vector summarizing the article counts per year and a a line plot of it.
#' @examples
#'\dontrun{
#' litReview(start = 2019, end = 2024, term = "Huntington's disease
#' AND histone post-translational modifications")
#'}
#'
#' @export

litReview <- function(start, end, term, db= "pubmed", save_plot = FALSE, save_file= FALSE){


  if (!requireNamespace("rentrez", quietly = TRUE)) {
    cli::cli_abort("Package 'rentrez' is required but not installed.")
  }

  # Years into which to look for publications
  year <- {{start}}:{{end}}


  # terms to search for in PubMed
  search_term <- stringr::str_glue("{term} AND {year}[PDAT]")


  df <- tidyr::tibble(year = year,
                      search = {{search_term}})


  #cat("\nCounting publications ...\n")

  df <- df |>
    dplyr::mutate(count = map_dbl({
      {
        search_term
      }
    }, ~ rentrez::entrez_search(
      db = {{db}},
      term = .x,
      use_history = TRUE
    )$count, .progress = TRUE))


  #cat("\nDone Counting publications.\n")


  #cat("\nScraping publication IDs ...\n")

  df <- df |>
    mutate(id= purrr::map({{search_term}}, ~rentrez::entrez_search(db = {{db}}, term= .x, use_history = TRUE)$id, .progress = TRUE)
    )
  #cat("\nDone scraping publication IDs.\n")


  total_queries <- df



  if(sum(df$count) == 0){

    cli::cli_alert_info('No publication was found!')
  }else{
  #cat("\nScraping publication titles ...\n")


  #function to get titles of articles




  title_df <- .get_titles(total_queries, db = {{db}})

  #cat("\nDone scraping publication titles!\n")


  start_year <- min(year)
  end_year <- max(year)


  title_df_long <- title_df |>
    dplyr::mutate(title = as.character(title)) |>
    dplyr::select(-c(search, count))


  title_df_wide <- title_df |>
    dplyr::mutate(title = as.character(title)) |>
    dplyr::group_by(year, search, count) |>
    dplyr::summarize(dplyr::across(id:title, ~paste0(na.omit(.x), collapse = ";")),
                     .groups = "drop")

  if(save_file == TRUE){

    write.csv(x = title_df_wide, str_glue("{term}_{start_year}_{end_year}.csv"), row.names = FALSE)
    write.csv(x = title_df_long, file = stringr::str_glue("long_{term}_{start_year}_{end_year}.csv"),
              row.names = FALSE)

    cat(stringr::str_glue("\n\nFile {term}_{start_year}_{end_year}.csv is saved successfully.\n"))

  }


  if(nrow(title_df)>0){

    df_plot <- title_df |> #dplyr::distinct() |>
      dplyr::select(Year = year, search) |>
      dplyr::count(Year, search, name = 'Count') |>
      dplyr::mutate(
                    search = stringr::str_remove_all(search, " AND \\d+\\[PDAT\\]"),
                    search = factor(search)
                    )

  p <- ggplot2::ggplot(df_plot, aes(x= Year, y= Count, color= search, group=1))+
      ggplot2::labs(color = "PubMed Search Term")+
      ggplot2::geom_line(linewidth= 1.2)+
      ggplot2::theme_classic(base_size = 12)+
      ggplot2::scale_color_brewer(palette = "Dark2")+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks())+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())+

    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size= 8)
    )
  if(save_plot){
    ggplot2::ggsave(filename = stringr::str_glue("{term}_{start_year}_{end_year}.png"),
                    width = 10, height = 7, dpi = 300, units = "in")
    }

  }else{p= NULL}
 return_list = list(Data = title_df_long, df_plot, Summary= table(title_df_long$year), plot= p)
 nonull_list = Filter(Negate(is.null), return_list)
return(nonull_list)


  }




}


#'@noRd
.get_titles <- function(df, db) {
df2 <- df |>
  dplyr::filter(count > 0) |>
  tidyr::unnest(id) |>
  dplyr::mutate(title = map_chr(
    .x = id,
    ~ rentrez::entrez_summary(db = db, id = .x)$title,
    .progress = TRUE
  ))

return(df2)

}
