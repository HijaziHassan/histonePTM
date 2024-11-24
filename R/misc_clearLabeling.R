



#' Reverse PTM Labeling
#'
#' For final representation it is better to replace chemical labeling with biological state. For example, replace 'propionyl' with 'unmodified'
#' or 'butyryl' with 'me1' (since me1+pr = bu).
#'
#' @param ptm_string PTM-containing string
#' @param rules A named vector containing \code{(labeled_ptm = unlabeled_ptm)}. Example: \code{c("pr" = "un")}.
#' @param labeling The labeling agent used. Use `PA` (propionic anhydride), `TMA` ( trimethylacetic anhydride), or `PIC_PA` (phenylisocynate + PA).
#' There are embedded rules for those reagents. If those don't match your case, define your own rules and `labeling`must not be `none`.
#' @param residue Specify wether to remove residue(capitlized letter followed by number). Either `keep` or `remove`.
#' @importFrom stringr str_detect str_replace_all
#' @importFrom cli cli_abort
#' @return unlabeled string or the same string if labeling is `none`.
#' @export
#'

misc_clearLabeling <- function(ptm_string, rules = NULL, labeling= NULL, residue = 'keep'){


  if(any(is.null(rules) && !labeling %in% c('none', 'PA', 'TMA', 'PIC_PA'))){
    cli::cli_abort('Please provide valid rules or specify a supported labeling argument.')
  }
  if(labeling!= "none"){

    if(is.null(rules)){

  nterm_rule = c(".+Nt-?" = "") #remove N-terminal Labeling

  prop_rules =  c(
    "pr"="un", #replace propionyl with unmod
    "bu"="me1" #replace bu which is me1-pr with me1
  )

  tma_rules = c(
    "tmame1" = "me1", #replace tmame1 with me1
    "tma" = "un") #replace TMAyl with unmod


  pic_rules =  c("pic" = "un") #replace phenylisocynate with unmod

  pic_pro_rules = c(prop_rules, pic_rules)

    }

  if(labeling== "PA"){

  if(any(stringr::str_detect(ptm_string, "(?:[:upper:]*\\d*me1)\\b")) & "me1" %in% unname(rules)){
    cli::cli_abort(c('You are trying to replace "bu" by "me1" but "me1" is already in your `rules`.',
                     "i" = "remove any 'me1' occurence in `ptm_string` before applying the function or change the `rules` argument."))
  }

 rules = c(nterm_rule, prop_rules)

  } else if (labeling== "TMA"){

    rules = c(nterm_rule, tma_rules)

  }else if(labeling== "PIC_PA"){

    rules = c(nterm_rule, prop_rules, pic_rules)
  }

    #if non of the labeling used, then the labeling arg is ignored and rules provided by the user is considered

unlabeld_string = stringr::str_replace_all(ptm_string, rules)

# option to keep or remove the residues which assumed to be capitalized and followed by a number
  if(residue == 'keep'){

  return(unlabeld_string)

} else if(residue == 'remove'){

stringr::str_remove_all(unlabeld_string, '[:upper:]\\d+')
    }
}else{return (ptm_string)}

}





