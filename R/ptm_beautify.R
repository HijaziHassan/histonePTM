#' ptm_beautify
#' @description It shorthands PTM names into pretty annotations
#' @param PTM The PTM string to be shorthanded
#' @param lookup The shorthanded names to replace the long PTM names (`histptm_lookup` and `shorthistptm_mass` are named vector provided by default.
#' The user can provide any named vector adapted to the PTM being used.)
#' @param software The software used for proteomic analysis (ex. Proline (Default) or Skyline)
#' @param residue Choose to keep or remove the residues. If `removed` ptms will be separated by a hyphen.
#'
#' @return renamed PTM string
#'
#'
#' @examples ptm_beautify(PTM = "TMAyl_correct (Any N-term); Butyryl (K28); Trimethyl (K37); Propionyl (K38)",
#'  lookup = histptm_lookup,
#'  software = "Proline", residue = 'keep')
#'
#'  ptm_beautify(PTM = "TMAyl_correct (Any N-term); Butyryl (K28); Trimethyl (K37); Propionyl (K38)",
#'  lookup = histptm_lookup,
#'  software = "Proline", residue = 'remove')
#'
#' ptm_beautify(PTM = "K[+126.06808]SAPATGGVK[+56.026215]K[+56.026215]PHR",
#' lookup = shorthistptm_mass,
#' software = "Skyline", residue = 'keep')
#'
#' @import stringr
#'


#function to round ptm masses. !!! this will remove (+) sign but keep (-) if exists





#' @export
ptm_beautify <- function(PTM,
                       lookup,
                       software = c("Proline", "Skyline"),
                       residue = c('keep', 'remove')){

  software <- match.arg(software)
  residue <- match.arg( residue)



  if (software == "Proline"){


      reversed_ptm <- stringr::str_replace_all(string = {{PTM}},
                                             c("([a-zA-Z0-9_\\:\\(\\)-]+\\s?[a-zA-Z0-9_\\:\\(\\)-]+?)\\s\\(([A-Z]\\d+)\\)"= "\\2\\1"))

      renamed_ptm <- stringr::str_replace_all(string = reversed_ptm, {{lookup}})

         if(residue == 'keep'){
             return(renamed_ptm)
      }  else if (residue == 'remove'){
      return(stringr::str_remove_all(renamed_ptm, '[:upper:]\\d+'))
      }


} else if (software == "Skyline"){

  # match what is inside the square brackets
  pattern <- "\\[(\\+\\d+\\.?\\d*)\\]"
  rounded_lookup <- setNames( {{lookup}}, sapply(names( {{lookup}}), .round_mod))

  if (residue == 'keep') {
  renamed_ptm <- stringr::str_replace_all({{PTM}}, pattern, function(mod) {
    # Extract what is inside the brackets
    ptm_mass <- stringr::str_match(mod, pattern)[2]
    ptm_mass_rounded <- as.character(.round_mod(ptm_mass))

    # Replace using the roundec lookup vector, if  exists
    if (ptm_mass_rounded %in% names(rounded_lookup) ) {
      return(rounded_lookup[[ptm_mass_rounded]])
    } else {
      return(mod)
    }
  })

  } else if (residue == 'remove'){


    # Find all matches
    matches <- str_match_all({{PTM}}, pattern)[[1]]

    # Extract numbers within brackets in column 2
    ptm_mass <- matches[, 2]

    #round the matched masses
    ptm_mass_rounded <- sapply(ptm_mass, .round_mod)

    # Replace the modifications using the normalized lookup vector
    replaced_ptm <- sapply(ptm_mass_rounded, function(mod) {
      mod_str <- as.character(mod)
      if (mod_str %in% names(rounded_lookup)) {
        return(rounded_lookup[[mod_str]])
      } else {
        return(mod)
      }
    })

    renamed_ptm <- paste(replaced_ptm, collapse = "-")


  }

return(renamed_ptm)

  }


}

#' @noRd
.round_mod <- function(mod) {
  round(as.numeric(mod), 1)
}

