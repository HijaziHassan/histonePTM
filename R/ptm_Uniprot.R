





#' Retreive MS-reported PTMs from Uniprot
#' @description
#' Get PTM information  reported in literature and documented by Uniprot.
#'
#'
#' @param Uniprot_accession Uniprot protein-specific accession code.
#' @param save_file \code{logical}. save results as a \code{.csv} file (optional).
#'
#' @import httr2
#' @import cli
#' @import tidyr
#'
#'
#' @return A \code{tibble} and \code{.csv} file.
#'
#'
#' @examples
#' ptm_Uniprot(Uniprot_accession = "B9FXV5", save_file = TRUE)
#'
#' @export

ptm_Uniprot <- function(Uniprot_accession, save_file= FALSE){

resp <- get_url(accession = {{Uniprot_accession}})


  resp_df <- resp |>
  httr2::resp_body_json() |>
   purrr::pluck("features") |>
  purrr::map_dfr(\(x) {as_tibble(x)})

  if('ptms' %in% names(resp_df)){

  resp_df <- resp_df |>
   tidyr::unnest_wider("ptms") |>
   tidyr::unnest_longer(c(dbReferences,sources)) |>
   tidyr::unnest_wider(c(xrefs, evidences, dbReferences), names_sep = "_") |>
   tidyr::unnest_wider(c(evidences_source, dbReferences_properties), names_sep = "_") |>
  conditional_unnest("evidences_source_properties")
  }else{

    cli::cli_alert_info("No PTMs are found for '{Uniprot_accession}' in Uniprot.")
}




  if(isTRUE(save_file) & length(resp_df > 0L)){
    file_name = paste0(Uniprot_accession, "_ptmUniprot.csv")
    write.csv(x = resp_df,file =  file_name, row.names = FALSE)
    cli::cli_alert_success("The file {file_name} is saved successfully.")
  }

  return(resp_df)
}




#' @noRd
get_url <- function(accession){
  PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"
  response_check <- httr2::request(PROTEINS_API) |>
    httr2::req_url_path_append('proteomics-ptm',
                               {{accession}}) |>
    httr2::req_error(is_error = \(resp) FALSE) |>
    httr2::req_perform()

  if(httr2::resp_status(response_check) == 400){
    cli::cli_abort(c("The request is failed.",
                     "x" = "The provided accession '{accession}'  is not valid.",
                     "i" = "Check the spelling or manually check it in Uniprot."))
  }

  return(response_check)
}



#' @noRd
conditional_unnest <- function(df, var, stopme= FALSE){

  ifelse({{var}} %in% names(df),
         return(tidyr::unnest_wider(df, {{var}})),
         return(df)
  )

}

