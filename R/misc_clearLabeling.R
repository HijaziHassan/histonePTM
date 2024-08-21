



#' Reverse PTM Labeling
#'
#' For final representation it is better to replace chemical labeling with biological state. For example, replace 'propionyl' with 'unmodified'
#' or 'butyryl' with 'me1' (since me1+pr = bu).
#'
#' @param ptm_string PTM-containing string
#' @param rules A named vector containing \code{(labeled_ptm = unlabeled_ptm)}. Example: \code{"pr" = "un"}.
#'
#' @return unlabeled string
#' @export
#'

misc_clearLabeling <- function(ptm_string, rules){

  if(is_missing(rules)){
  rules = c(".+Nt-?" = "", #remove N-terminal Labeling
            "pr"="un", #replace propionyl with unmod
            "bu"="me1", #replace bu which is me1-pr with me1
            "tmame1" = "me1", #replace tmame1 with me1
            "tma" = "un" #replace TMAyl with unmod
  )

  if(any(stringr::str_detect(ptm_string, "(?:[:upper:]*\\d*me1)\\b")) & "me1" %in% unname(rules)){
    cli::cli_abort(c('You are trying to replace "bu" by "me1" but "me1" is already in your `rules`.',
                     "i" = "remove any 'me1' occurnce in `ptm_string` before applying the function or change the `rules` argument."))
  }

  }


    stringr::str_replace_all(ptm_string, rules)

}





