


#extract ptm and its index
.extract_ptm_indx <- function(mod){

  stringr::str_match_all(string = mod,
                         pattern = "(?<modif>[a-zA-Z0-9_-]+\\s?[a-zA-Z0-9_-]+?) \\([A-Z](?<index>\\d+)\\)") |>
    data.frame() |>
    dplyr::mutate(modif = paste0("[",stringr::str_replace_all(modif, histptm_mass), "]"))
}

#insert ptm in its respective position in the sequence
.insert_ptm_seq <- function(seq, replace){

  if (nrow(replace) == 0) { #if no modification, return the bare sequence
    return(seq)
  }

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
  return(seq)
}




#' Proforma modified peptide representation
#'
#' @param seq The sequence to be modified.
#' @param mod The modifications with their positions on the peptide
#'
#' @return sequence with monoisotopic mass of ptm at its modification site.
#'
#' @import purrr
#' @import stringr
#' @import purrr
#' @examples
#' ptm_toProforma(seq ="KSAPATGGVKKPHR", mod= "Propionyl (Any N-term); Lactyl (K1); Dimethyl (K10); Propionyl (K11)")
#'
#' @examples
#' ptm_toProforma(seq ="KQLATKVAR", mod= "Propionyl (Any N-term); Propionyl (K1); Propionyl (K6)", Nterm = "[+56.026]")
#'
#' @export
ptm_toProforma <- function(seq, mod, Nterm = FALSE){

modified_peptide <- purrr::map2_chr(

    .x = seq,
    .y =  purrr::map(mod, .extract_ptm_indx),
    .f =~ .insert_ptm_seq(seq = .x, replace = .y),
    .progress = TRUE)
if(Nterm == FALSE){
return(noquote(modified_peptide))}else{
return(noquote(paste0(Nterm,"-", modified_peptide)))
}

}




