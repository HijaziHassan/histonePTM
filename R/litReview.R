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
#' @param save_file \code{logical}. If `TRUE` results will be saved as `.csv` files.
#' @import rentrez
#' @importFrom stringr str_glue
#' @importFrom purrr map
#' @importFrom dplyr mutate select summarize filter group_by across
#' @importFrom tidyr tibble unnest
#' @importFrom utils write.csv
#'
#' @return A list. A dataframe with  year, id, and article title. A vector summarizing the article counts per year.
#' @examples
#'\dontrun{
#' litReview(start = 2019, end = 2024, term = "Huntington's disease
#' AND histone post-translational modifications")
#'}
#'
#' @export

litReview <- function(start, end, term, db= "pubmed", save_file= FALSE){

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

return(list((title_df_long), table(title_df_long$year)))


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
