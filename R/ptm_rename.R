#' ptm_rename
#' @description It shorthands PTM names into pretty annotations
#' @param PTM The PTM string to be shorthanded
#' @param lookup The shorthanded names to replace the long PTM names (histptm_lookup is a named vector provided by default. The user can provide any named vector adapted to the PTM being used.)
#' @param software The software used for proteomic analysis (ex. Proline (Default), MaxQaunt, Skyline, Fragpipe, ...)
#'
#' @return renamed PTM string
#'
#'
#' @examples ptm_rename(PTM = "TMAyl_correct (Any N-term); Butyryl (K28); Trimethyl (K37); Propionyl (K38)",
#'  lookup = histptm_lookup,
#'  software = "Proline")
#'
#'
#' @import stringr
#' @export ptm_rename
#'
ptm_rename <- function(PTM, lookup = histonePTM::histptm_lookup, software = "Proline"){

  if (software == "Proline"){

    reversed_ptm <- stringr::str_replace_all(string = {{PTM}}, c("([a-zA-Z0-9_\\:\\(\\)-]+\\s?[a-zA-Z0-9_\\:\\(\\)-]+?)\\s\\(([A-Z]\\d+)\\)" = "\\2\\1"))

    renamed_ptm <- stringr::str_replace_all(string = reversed_ptm, {{lookup}})

    return(renamed_ptm)

  }
}


