#' Shift position numbers in histone modification labels
#'
#' In histone literature, amino acid positions are conventionally numbered starting
#' from the first residue after the cleaved N-terminal methionine (iMet).
#' However, in databases like UniProt, the initiator methionine is often included
#' in the sequence, resulting in a +1 offset. This function shifts the numbers
#' to convert between these numbering systems (e.g., converting UniProt positions
#' to the histone field standard).
#'
#' By default, it shifts numbers immediately preceded by **K**, **S**, or **T**
#' (the most common modified residues in histones). You can supply a custom
#' pattern for other residues (e.g., R).
#'
#' @param labels Character vector of labels containing numbers to shift.
#' @param by Integer value to add to each number (default: -1).
#'        Use -1 to convert from UniProt numbering (including Met) to the
#'        histone field standard (excluding Met). Use +1 for the reverse.
#' @param pattern Regular expression pattern for the numbers to shift.
#'        Default is \code{"(?<=[KST])[0-9]+"} (numbers immediately preceded
#'        by K, S, or T).
#'
#' @return Character vector with shifted numbers.
#' @export
#'
#' @examples
#' # Shift lysine positions (remove Met offset)
#' misc_shiftNumbers("tmaNt-K19ac-K24tma", by = -1)
#' # Result: "tmaNt-K18ac-K23tma"
#'
#' # Shift serine/threonine positions (phosphorylation sites)
#' misc_shiftNumbers("S10ph-T3ph", by = -1)
#' # Result: "S9ph-T2ph"
#'
#' # Works with semicolon‑separated lists too
#' misc_shiftNumbers("Propionyl (Any N-term); Dimethyl (K10); Propionyl (K15)", by = -1)
#' # Result: "Propionyl (Any N-term); Dimethyl (K9); Propionyl (K14)"
#'
#' # Convert from histone standard back to UniProt numbering (add Met offset)
#' misc_shiftNumbers(c("K27me3", "S10ph"), by = 1)
#'
#' # Custom pattern: shift numbers after R (arginine)
#' misc_shiftNumbers("R2me1", pattern = "(?<=R)[0-9]+", by = -1)
misc_shiftNumbers <- function(labels, by = -1, pattern = "(?<=[KST])[0-9]+") {
  stringr::str_replace_all(labels, pattern, function(m) {
    nums <- as.numeric(m)
    # nums are always numeric;
    # but this is a safe fallback (no warning)
    ifelse(is.na(nums), m, as.character(nums + by))
  })
}
