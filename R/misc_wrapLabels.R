#' Wrap modification labels for plotting
#'
#' Wraps histone PTM labels so they fit within plot margins.
#' Handles three common formats:
#' \itemize{
#'   \item Semicolon-separated lists (e.g., "K4me1; K9me3") — splits at "; " and
#'         attempts to keep each segment under \code{max_width}.
#'   \item Hyphen/minus-chained labels (e.g., "K27me3-K36me2") — wraps at hyphens,
#'         keeping the hyphen attached to the preceding chunk.
#'   \item Plain text — simple word wrapping on spaces.
#' }
#'
#' @param ptm_strings Character vector of PTM labels to wrap.
#' @param max_width Maximum number of characters per line before wrapping.
#'        Default is 20.
#' @importFrom purrr reduce map_chr
#' @importFrom stringr str_detect str_split
#' @return A character vector of the same length as \code{ptm_strings}, with
#'         newline characters (\code{\\n}) inserted where wrapping occurs.
#'
#' @export
#'
#' @examples
#' # Semicolon-separated list
#' misc_wrapLabels("Propionyl (Any N-term); Dimethyl (K10); Propionyl (K15)")
#'
#' # Hyphen-chained labels
#' misc_wrapLabels("tmaNt-K19ac-K24tma")
#'
#' # Force narrower wrapping
#' misc_wrapLabels("tmaNt-K19ac-K24tma", max_width = 10)
#'
#' # Mixed vector
#' labels <- c(
#'   "tmaNt-K19ac-K24tma",
#'   "Propionyl (Any N-term); Dimethyl (K10); Propionyl (K15)"
#' )
#' cat(misc_wrapLabels(labels), sep = "\n\n")
misc_wrapLabels <- function(ptm_strings, max_width = 20) {

  # Helper for hyphen/minus chains
  wrap_hyphen_chain <- function(hyphenated_ptm_string, max_width) {
    chunks <- stringr::str_split(hyphenated_ptm_string, "(?<=[-\u2212])")[[1]]
    chunks <- chunks[nchar(chunks) > 0]

    if (length(chunks) == 0) return("")

    state <- purrr::reduce(
      .x = chunks,
      .f = function(acc, chunk) {
        if (acc$current == "") {
          return(list(lines = acc$lines, current = chunk))
        }
        candidate <- paste0(acc$current, chunk)
        if (nchar(candidate) > max_width) {
          list(lines = c(acc$lines, acc$current), current = chunk)
        } else {
          list(lines = acc$lines, current = candidate)
        }
      },
      .init = list(lines = character(0), current = "")
    )

    all_lines <- c(state$lines, state$current)
    all_lines <- all_lines[nchar(all_lines) > 0]
    paste(all_lines, collapse = "\n")
  }

  # Helper for semicolon lists
  wrap_semicolon <- function(label, max_width) {
    segs <- stringr::str_split(label, "; ")[[1]]
    segs <- segs[nchar(segs) > 0]

    if (length(segs) == 0) {
      return(label)
    }

    state <- purrr::reduce(
      .x = seq_along(segs)[-1],
      .f = function(acc, i) {
        next_seg <- segs[i]
        candidate <- paste0(acc$current, "; ", next_seg)
        if (nchar(candidate) > max_width) {
          new_line <- paste0(acc$current, ";")
          list(lines = c(acc$lines, new_line), current = next_seg)
        } else {
          list(lines = acc$lines, current = candidate)
        }
      },
      .init = list(lines = character(0), current = segs[1])
    )

    all_lines <- c(state$lines, state$current)
    all_lines <- all_lines[nchar(all_lines) > 0]
    paste(all_lines, collapse = "\n")
  }

  # Main vectorized loop
  purrr::map_chr(
    .x = ptm_strings,
    .f = function(ptm_string) {
      if (is.na(ptm_string)) return(NA_character_)
      if (stringr::str_detect(ptm_string, "; ")) {
        wrap_semicolon(ptm_string, max_width)
      } else if (stringr::str_detect(ptm_string, "[-\u2212]")) {
        wrap_hyphen_chain(ptm_string, max_width)
      } else {
        # Fallback: hard-max wrapper on spaces
        words <- stringr::str_split(ptm_string, "\\s+")[[1]]
        state <- purrr::reduce(
          .x = words,
          .f = function(acc, w) {
            if (acc$current == "") {
              return(list(lines = acc$lines, current = w))
            }
            candidate <- paste0(acc$current, " ", w)
            if (nchar(candidate) > max_width) {
              list(lines = c(acc$lines, acc$current), current = w)
            } else {
              list(lines = acc$lines, current = candidate)
            }
          },
          .init = list(lines = character(0), current = "")
        )
        all_lines <- c(state$lines, state$current)
        all_lines <- all_lines[nchar(all_lines) > 0]
        paste(all_lines, collapse = "\n")
      }
    }
  )
}
