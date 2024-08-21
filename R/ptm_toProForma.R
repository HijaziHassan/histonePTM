


#extract ptm and its index
.extract_ptm_indx <- function(mod){

  stringr::str_match_all(string = mod,
                         pattern = "(?<modif>[a-zA-Z0-9_-]+\\s?[a-zA-Z0-9_-]+?) \\([A-Z](?<index>\\d+)\\)") |>
    data.frame() |>
    dplyr::mutate(modif = paste0("[", stringr::str_replace_all(modif, histonePTM::histptm_mass), "]"))
}

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




#' Convert sequence-PTM to ProForma notation.
#' @description
#' Convert Proline modification string to ProForma notation.
#'
#'
#' @param seq The sequence to be modified.
#' @param mod The modifications with their positions on the *peptide*.
#' @param Nterm N-terminal modification
#'
#' @return character string
#'
#' @importFrom purrr map map2_chr
#' @importFrom stringr str_detect str_match_all str_replace_all
#' @importFrom dplyr mutate
#' @examples
#' ptm_toProForma(seq = "KSAPATGGVKKPHR",
#'                mod = "Propionyl (Any N-term); Lactyl (K1); Dimethyl (K10); Propionyl (K11)")
#'
#' @examples
#' ptm_toProForma(  seq = "KQLATKVAR",
#'                  mod = "Propionyl (Any N-term); Propionyl (K1); Propionyl (K6)",
#'                Nterm = "[+56.026]")
#'
#' @export

ptm_toProForma <- function(seq, mod, Nterm = FALSE){

if(any(stringr::str_detect(mod, pattern = 'Propionyl \\(Any N-term\\)'))){
  Nterm = "[+56.026]"
}else if(any(stringr::str_detect(mod, pattern = 'TMAyl_correct \\(Any N-term\\)'))){
  Nterm = "[+84.057515]"
}else{Nterm = Nterm}

modified_peptide <- purrr::map2_chr(

    .x = seq,
    .y = purrr::map(mod, .extract_ptm_indx),
    .f = ~ .insert_ptm_seq(seq = .x, replace = .y),
    .progress = TRUE)


if(Nterm == FALSE){
return(noquote(modified_peptide))}else{
return(noquote(paste0(Nterm,"-", modified_peptide)))
}

}




