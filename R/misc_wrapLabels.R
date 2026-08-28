#' Wrap modification labels for plotting
#'
#' This function wraps histone modification labels to fit within plot margins.
#' It handles three common formats:
#' \itemize{
#'   \item Semicolon-separated lists (e.g., "K4me1; K9me3") → split into separate lines
#'   \item Hyphen/minus-chained labels (e.g., "K27me3-K36me2-K79me1") → wrap at hyphens,
#'         keeping the hyphen attached to the preceding chunk
#'   \item Plain text → standard word wrapping using \code{stringr::str_wrap()}
#' }
#'
#' @param labels Character vector of labels to wrap.
#' @param max_width Maximum number of characters per line before wrapping.
#'        Default is 20.
#'
#' @return A character vector of the same length as \code{labels}, with newline
#'         characters (\code{"\\n"}) inserted where wrapping occurs.
#'
#' @importFrom stringr str_split str_detect str_replace_all str_wrap
#' @importFrom purrr reduce map_chr
#'
#' @export
#'
#' @examples
#' # Example labels from real proteomics data
#' labels <- c(
#'   "TMAyl_correct (Any N-term); TMAyl-me1_DP (K1); Acetyl (K6)",
#'   "tmaNt-K28tmame1-K37tma-K38tma",
#'   "tmaNt-me3-tmame1-tma"
#' )
#'
#' # Apply default wrapping (max_width = 20)
#' wrapped <- misc_wrapLabels(labels)
#'
#' # View raw strings with \n characters
#' print(wrapped)
#'
#' # View formatted output with actual line breaks
#' cat("Formatted labels:\n")
#' cat(wrapped, sep = "\n\n")
#'
#' # Try a narrower width to force more breaks
#' wrapped_narrow <- misc_wrapLabels(labels, max_width = 10)
#' cat("\nNarrow width (10 chars):\n")
#' cat(wrapped_narrow, sep = "\n\n")
#'
#' # Works with single labels too
#' misc_wrapLabels("K27me3-K36me2-K79me1-K80me2", max_width = 12)
#'
misc_wrapLabels <- function(labels, max_width = 20) {

  wrap_hyphen_chain <- function(label, max_width) {
    Kptm <- stringr::str_split(label, "(?<=[-\u2212])")[[1]]

    final_state <- purrr::reduce(
      .x = Kptm,
      .f = function(acc, p) {
        candidate <- paste0(acc$current, p)
        if (nchar(candidate) > max_width && acc$current != "") {
          list(lines = c(acc$lines, acc$current), current = p)
        } else {
          list(lines = acc$lines, current = candidate)
        }
      },
      .init = list(lines = character(0), current = "")
    )

    all_lines <- c(final_state$lines, final_state$current)
    all_lines <- all_lines[nchar(all_lines) > 0]
    paste(all_lines, collapse = "\n")
  }


  purrr::map_chr(
    .x = labels,
    .f = function(label) {
      if (stringr::str_detect(label, "; ")) {
        stringr::str_replace_all(label, "; ", "\n")
      } else if (stringr::str_detect(label, "[-\u2212]")) {
        wrap_hyphen_chain(label, max_width)
      } else {
        stringr::str_wrap(label, width = max_width)
      }
    }
  )
}
