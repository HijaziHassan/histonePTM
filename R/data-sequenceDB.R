#' @title Histone protein and peptide sequence labeling
#' @description
#' a named list of proteins including a nested list of sequence, enzyme fragments and labels
#'
#' \itemize{
#'   \item Protein. H3, H4, H2A, H2B and iRT
#'      \item seq. Full amino acid sequence of the core protein and its variants.
#'      \item Fragments. argC_frags/tryp_frags.
#'      \item label. Representation of fragments by position in sequence (e.g. "H3 (3-8)" )
#'
#' }
#' @format a named character vector of proteins.
#' @source in-house user defined names
#' @examples
#' # sequenceDB[['H3']]
#' # sequenceDB[['H4']][['label']]
#' # sequenceDB$H3$argC_frags
#' # sequenceDB$H3$seq
#'

#'
#'
#'
#'

"sequenceDB"
