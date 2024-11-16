


#' Convert sequence-PTM to ProForma notation.
#' @description
#' Convert modified peptide sequence to ProForma format.
#'
#'
#' @param seq The stripped sequence to be modified (if modification are already embeded in the sequence, set `replace_only` to `TRUE` and no need to provide `mod`.
#' @param mod The modifications with their positions on the *peptide* (as represented in `Proline` software)
#' @param lookup A named vector contaning PTMs and their replacment values (such as `histptm_mass` or `histptm_unimod` (`default`) or any custom vector). DON'T use names as "GG" as modifications. This will replace also "GG" residues if found in the modified sequence provided in `mod`.
#' @param replace_only (\code{Optional}). If `TRUE`, it assumes that no need for a `seq` argument. It only replaces values in \code{mod} according to the provided `lookup` vector.
#' @param Nterm N-terminal modification. Add systematically custom "Nterm" modification to all the sequences.
#'
#' @return character string
#'
#' @importFrom purrr map map2_chr
#' @importFrom stringr str_detect str_match_all str_replace_all str_extract
#' @importFrom dplyr mutate
#' @examples
#' ptm_toProForma(seq = "KSAPATGGVKKPHR",
#'                mod = "Propionyl (Any N-term); Lactyl (K1); Dimethyl (K10); Propionyl (K11)"
#'                )
#'
#' @examples
#' ptm_toProForma(seq = "KQLATKVAR",
#'                mod = "Propionyl (Any N-term); Propionyl (K1); Propionyl (K6)",
#'                lookup = histptm_mass
#'                )
#' @examples
#' ptm_toProForma(seq = 'KacSAPATGGVKprKprPHR',
#'               lookup = c(ac = 'UNIMOD:1', pr= 'UNIMOD:58'),
#'               replace_only = TRUE
#'               )
#'
#' @export

ptm_toProForma <- function(seq, mod, lookup = NULL, replace_only = FALSE, Nterm = ""){


if(replace_only == TRUE){

  if(is.null(lookup)) lookup = setNames(names(shorthistptm_mass), shorthistptm_mass)
  #add square brackets for any vector used for replacement  if the values are not already in between []
  if(any(grepl('\\[', lookup)) == FALSE){

  lookup = sapply(X = lookup, \(X) paste0("[", X, "]"))}

  modified_peptide =  stringr::str_replace_all(seq, lookup)

}else{

  if(is.null(lookup)) lookup = histptm_unimod

  modified_peptide <- purrr::map2_chr(

    .x = seq,
    .y = purrr::map(.x = mod, .f = ~.extract_ptm_indx(.x, lookup = lookup)),
    .f = ~ .insert_ptm_seq(seq = .x, replace = .y),
    .progress = TRUE)

  nterm <- .extract_nterm(mod)

  if(!is.na(nterm)){

    nterm = stringr::str_replace_all(nterm, lookup)

    modified_peptide <- paste0("[", nterm, "]-", modified_peptide)
  }else if(Nterm != ""){
    modified_peptide <- paste0(Nterm,"-", modified_peptide)
  }

      }

  modified_peptide

}


#' @noRd
#extract ptm and its index
.extract_ptm_indx <- function(mod, lookup){

  stringr::str_match_all(string = mod,
                         pattern = "(?<modif>[a-zA-Z0-9_-]+\\s?[a-zA-Z0-9_-]+?) \\([A-Z](?<index>\\d+)\\)") |>
    data.frame() |>
    dplyr::mutate(modif = paste0("[", stringr::str_replace_all(modif, lookup), "]"))
}

#' @noRd
#insert ptm in its respective position in the sequence
.insert_ptm_seq <- function(seq, replace){

  if (nrow(replace) != 0) { #if no modification, return the bare sequence


    replacement <- as.character(replace$modif)
    position <- as.integer(replace$index)

    #move backward
    for(i in nrow(replace):1){
      upstream <- substr(seq,
                         start = 1L,
                         stop = position[[i]])

      downstream <- substr(seq,
                           start = position[[i]] + 1L,
                           stop = nchar(seq))

      # update sequence
      seq <- paste0(upstream, replacement[[i]], downstream)
    }

  }
  return(seq)
}

#' @noRd
.extract_nterm <- function(mod){

  stringr::str_extract(string = mod,
                       pattern = '(.+)\\s\\(Any N-term\\)', group = 1)
}

