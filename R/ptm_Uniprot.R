
# Helper function to download data
get_url <- function(url, ...) {
  response <- httr::GET(url, ...)

  if (httr::http_status(response)$category != "Success") {
    cat(httr::content(response, as = "text"))
    stop("Error in HTTP request:", httr::http_status(response)$reason)
  }

  return(response)
}


PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"


#' Retreive proteomics-validated PTMs from Uniprot
#' @description
#' A function that access to Uniprot to retreive documented and classifed PTMs along with
#' many other features such as universal spectrum idetifier (USI) if reported by the database of origin.
#'
#' @param accession The accession number of the interrogated protein.
#' @param config For advanced users.
#' @param ... for any later addition.
#' @import rlang
#' @import stringr
#' @import tidyr
#' @import httr
#' @import jsonlite
#'
#' @return A \code{tibble} with PTMs reported in Uniprot with their labels (gold, silver, ...).
#' Note not all PTMs are yet documented in Uniport.
#' @examples
#' ptm_Uniprot(acession = "B9FXV5")
#' @export

ptm_Uniprot <- function(accession, config = httr::accept("application/json"), ...){

  accession = rlang::enexpr(accession)

  full_url= stringr::str_glue("{PROTEINS_API}/proteomics-ptm/{accession}")

  #select the columns that I think would be of interest.
  # cols_to_keep <- c('accession', 'entryName', 'taxid', 'begin', 'end', 'peptide',
  #                   'name', 'position', 'sources', 'id', 'properties.Pubmed ID',
  #                   'properties.Confidence score','properties.Universal Spectrum Id',
  #                   'properties.Proforma', 'dbReferences.properties.Localization probability')


  response = get_url(url= {{full_url}})

  #turn json to a readable format
  data_chr <- base::rawToChar(response$content)

  #then to a dataframe
  data_json <-jsonlite::fromJSON(data_chr, flatten = TRUE)

  #successive unnesting to only extract data of interest
  df_final <- data_json |>
    tidyr::as_tibble() |>
    #tidyr::unnest_auto(c(features, ptms, dbReferences)) |>
    #tidyr::unnest(cols = c(features, ptms, dbReferences)) |>
    tidyr::unnest(features) |>
    tidyr::unnest(ptms) |>
    tidyr::unnest_longer(dbReferences) |>
    tidyr::unpack(dbReferences, names_sep = ".")
    #tidyr::unnest_wider(dbReferences)
    #tidyr::unnest_longer(dbReferences)
    # tidyr::unnest(dbReferences,
    #               names_sep = "." ) |>
    #  dplyr::select(any_of(cols_to_keep), ...)

  return(df_final)
}





