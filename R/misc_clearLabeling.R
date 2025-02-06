



#' Reverse PTM Labeling
#'
#' For final representation it is better to replace chemical labeling with biological state. For example, replace 'propionyl' with 'unmodified'
#' or 'butyryl' with 'me1' (since me1+pr = bu).
#'
#' @param ptm_string PTM-containing string
#' @param rules A named vector containing \code{(labeled_ptm = unlabeled_ptm)}. Example: \code{c("pr" = "un")}.
#' @param labeling The labeling agent used. Use `PA` (propionic anhydride), `TMA` ( trimethylacetic anhydride), or `PIC_PA` (phenylisocynate + PA).
#' If none applies or you don't want to use the built-in rules, set it to 'none' and indicate your `rules`. Otherwise, it defaults to `NULL`, returning the input as is.
#' There are embedded rules for those reagents. If those don't match your case, define your own rules and `labeling`must not be `none`.
#' @param residue Specify wether to remove residue(capitlized letter followed by number). Either `keep` or `remove`.
#' @importFrom stringr str_detect str_replace_all
#' @importFrom cli cli_abort cli_alert_warning
#' @return unlabeled string or the same string if labeling is `none`.
#' @export
#'

misc_clearLabeling <- function(ptm_string, rules = NULL, labeling = NULL, residue = c("keep", "remove")) {

  residue <- match.arg(residue)

  # Check if both rules and labeling are NULL
  if (is.null(rules) && is.null(labeling)) {
    cli::cli_alert_warning('Please provide valid rules or specify a supported labeling argument.')
    return(ptm_string)
  }


  supported_labels <- c('none', 'PA', 'TMA', 'PIC_PA')

  # Validate labeling argument
  if (!is.null(labeling)) {
    if (!labeling %in% supported_labels) {
      cli::cli_abort(c('x' = 'The provided labeling "{labeling}" is invalid.',
                       'i' = 'Please choose one of the following: {toString(supported_labels)}.'))
    }
  }

  # Default rules for different labeling types
  nterm_rule <- c(".+Nt-?" = "")  # Remove N-terminal labeling
  prop_rules <- c(
    "pr" = "un",  # Replace propionyl with unmodified
    "bu" = "me1"  # Replace butyrylation (me1-pr) with me1
  )
  tma_rules <- c(
    "tmame1" = "me1",  # Replace tmame1 with me1
    "tma" = "un"       # Replace TMAyl with unmodified
  )
  pic_rules <- c("pic" = "un")  # Replace phenylisocyanate with unmodified

  # Combine rules based on labeling
  if (!is.null(labeling)) {
    if (is.null(rules)) {
      if (labeling == "PA") {
        # Validate conflicting rules
        if (any(stringr::str_detect(ptm_string, "\\b(me1|[:upper:]\\d*me1)\\b")) && "me1" %in% unname(prop_rules)) {
          # cli::cli_abort(c(
          #   'x' = 'You are trying to replace "bu" with "me1", but "me1" is already in your rules.',
          #   'i' = 'Remove any "me1" occurrences in `ptm_string` before applying the function, or adjust the `rules` argument.'
          # ))
          cli::cli_alert_danger("'me1' exist in your PTMs and you are converting 'bu' to 'me1'.
                                'bu' will not be converted to 'me1'.")
          prop_rules <- prop_rules[-2]
        }

        rules <- c(nterm_rule, prop_rules)

      } else if (labeling == "TMA") {
        if (any(stringr::str_detect(ptm_string, "\\b(me1|[:upper:]\\d*me1)\\b")) && "me1" %in% unname(prop_rules)) {
          cli::cli_alert_danger("'me1' exist in your PTMs and you are converting 'tmame1' to 'me1'.
                                'tmame1' will not be converted to me1.")
          tma_rules <- tma_rules[-1]
        }

        rules <- c(nterm_rule, tma_rules)

      } else if (labeling == "PIC_PA") {
        rules <- c(nterm_rule, prop_rules, pic_rules)
      } else if (labeling == "none") {
        # Ignore labeling and use user-provided rules
        unlabeld_string <- stringr::str_replace_all(ptm_string, rules)
      }
    }
  }

  # Apply the rules to the ptm_string
  if (!is.null(rules)) {
    unlabeld_string <- stringr::str_replace_all(ptm_string, rules)
  } else {
    unlabeld_string <- ptm_string
  }

  # remove Capital letter(i.e. residue) followed by a number (i.e its position).
  if (residue == 'remove') {
    return(stringr::str_remove_all(unlabeld_string, '[:upper:]\\d+'))
  }

  # Return the unlabelled string
  return(unlabeld_string)
}





