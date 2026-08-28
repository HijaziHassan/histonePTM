#' Shift position numbers in histone modification labels
#'
#' In histone literature, amino acid positions are conventionally numbered starting
#' from the first residue after the cleaved N-terminal methionine (iMet).
#' However, in databases like UniProt, the initiator methionine is often included
#' in the sequence, resulting in a +1 offset. This function shifts the numbers
#' to convert between these numbering systems (e.g., converting UniProt positions
#' to the histone field standard).
#'
#' @param labels Character vector of labels containing numbers to shift.
#' @param by Integer value to add to each number (default: -1).
#'        Use -1 to convert from UniProt numbering (including Met) to the
#'        histone field standard (excluding Met). Use +1 for the reverse.
#' @param pattern Regular expression pattern for the numbers to shift.
#'        Default is \code{"(?<=K)[0-9]+"} (numbers immediately preceded by "K").
#'
#' @return Character vector with shifted numbers.
#' @export
#'
#' @examples
#' # Convert from UniProt numbering to histone standard (remove Met offset)
#' misc_shiftNumbers("tmaNt-K19ac-K24tma", by = -1)
#' # Result: "tmaNt-K18ac-K23tma"
#'
#' # Works with semicolon‑separated lists too
#' misc_shiftNumbers("Propionyl (Any N-term); Dimethyl (K10); Propionyl (K15)", by = -1)
#' # Result: "Propionyl (Any N-term); Dimethyl (K9); Propionyl (K14)"
#'
#' # Convert from histone standard back to UniProt numbering (add Met offset)
#' misc_shiftNumbers(c("K27me3", "K36me2"), by = 1)

misc_shiftNumbers <- function(labels, by = -1, pattern = "(?<=K)[0-9]+") {
  stringr::str_replace_all(labels, pattern, function(m) {
    num <- as.numeric(m)
    if (is.na(num)) {
      warning("Non-numeric match found: ", m, " — returning as-is")
      return(m)
    }
    as.character(num + by)
  })
}
