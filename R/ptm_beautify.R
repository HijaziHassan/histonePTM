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
  if(rlang::is_missing(lookup)){
    cli::cli_abort(c('You forget to add "lookup" argument.',
      'i' = 'provide a named vector depending on the PTMs you want to repalce.
      Ex: c("oldptm" = "new ptm").',
      'x'= '"lookup" argument is missing.'))
  }


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

  rounded_lookup <- setNames( {{lookup}}, sapply(names( {{lookup}}), .round_mod))

  if (residue == 'keep') {
    pattern <- "\\[(\\+\\d+\\.?\\d*)\\]"

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

#if the Skyline mdofied string was the one of the three letter code [2Me] for me2 ...
#then it returns the same string but with the square brackets.
if(renamed_ptm == {{PTM}}){
  renamed_ptm <- str_remove_all({{PTM}}, '\\[|\\]')
}



  } else if (residue == 'remove'){

    #pattern <- '[A-Z]?(?<=\\[)([-?+?0-9a-zA-Z]+)(?=\\])[A-Z]?'
    pattern <- '[A-Z]?(?<=\\[)([-+?0-9a-zA-Z]*\\.?[0-9a-zA-Z]+)(?=\\])[A-Z]?'

    renamed_ptm <-  purrr::map_chr(stringr::str_extract_all({{PTM}}, pattern = pattern),
                   ~ stringr::str_c(.round_mod(.x), collapse="-")
                   ) |>
      stringr::str_replace_all(rounded_lookup)



  }

  #Sklyline doesn't accept two modifications on the same residue.
  # for Nterm+Kmod we make custome modifications like [+72.021129] which represents la-Nt.
  #When replacing the string of sequence it nice to have Nt- before the sequence not after N-term K.
  #This if statement removes it infort of N-terminal K and return back.
  if(stringr::str_detect(string = renamed_ptm, pattern = 'prNt-|tmaNt-')){
    Nterm <- stringr::str_extract(string = renamed_ptm, pattern = 'prNt-|tmaNt-')
    renamed_ptm <- stringr::str_remove(string = renamed_ptm, pattern = Nterm)
    renamed_ptm <- paste0(Nterm, renamed_ptm)}

  return(renamed_ptm)


  }


}

#' @noRd
.round_mod <- function(mod) {
  if(anyNA(suppressWarnings(as.numeric(mod)))){
    return(mod)
  }else{
  round(as.numeric(mod), 1)}
}

