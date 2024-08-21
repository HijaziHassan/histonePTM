



#flag redundant isoforms having the same peak intensities.
#' Flagging suspicious sequence-PTM combinations
#' @description
#' In Proline result output, some sequence-PTM combinations have identical intensity values in all samples.
#' Usually those are isobaric PTM-sequence combinations that are most likely false positive hits.
#'
#'
#' @param df A \code{dataframe}.
#' @param PTM_col PTM column name in the dataset.
#' @param var_col Column name of a unique identifier (e.g. peptide sequence) in the dataset.
#' @param select_cols Either one of the select helper functions with proper argument
#' (e.g. \code{starts_with("WT_")}) or a vector \code{c()} containing the names of
#'  intensity/abundance columns to check for duplication.
#'
#' @import stringr
#' @import dplyr
#' @import rlang
#' @import cli
#'
#' @return A named \code{list} of two \code{tibbles}.
#'
#' \emph{unique} \cr
#' * a non-redundant version of \code{df} after:\cr
#'     - removing identical duplicate rows \cr
#'     - removing redundant rows sharing identical intensity values \cr
#'     - marking with asterisk (*) redundant PTMs in \code{PTM_col}.\cr
#'       Redundant observations are reported in \code{duplicate} (see below) with their frequency \code{n}.
#'
#' \emph{duplicate} \cr
#' * A subset \code{tibble} of \code{df} containing all de-duplicated identifications sharing
#'  identical intensities in \strong{ALL} \code{select_cols} columns.
#'
#' \emph{freq_table} \cr
#' * A table summarizing the duplicate counts and their frequencies in the dataset.
#'
#' @export

ptm_flagDupes <- function(df,PTM_col, var_col, select_cols) {

  select_exp <- rlang::enexpr(select_cols)
  PTM_col <-  rlang::enexpr(PTM_col)
  var_col <-  rlang::enexpr(var_col)



  # Record time
  start_time <- Sys.time()

  #remove IDENTICAL ROWS## could be same peptide from different proteins
  NotIdentical <- df[!duplicated(df), ]
  ##count the number of duplicate row(s) deleted##
  nrows_ident <- df[duplicated(df), ] |> nrow()
  cli::cli_alert_info('Identical rows removed: {nrows_ident}.')



  #isolate different IDs but having same intensities over ALL samples.
  Duplicates <- NotIdentical |>
    dplyr::add_count(dplyr::across(!!select_exp)) |>
    dplyr::filter(n > 1)

  rlang::eval_tidy(Duplicates)


  nrows_dup <- nrow(Duplicates)

  cli::cli_alert_info('Rows with duplicate intensities: {nrows_dup}.')


  #get the frequence of duplicate counts
freq_table <- table(Duplicates$n) |> as.data.frame()

if(nrow(freq_table)>0){
names(freq_table) <- c('dupes_count', 'freq')
}

  #remove duplications where:
  # the intensity value is identical for more than one ID (case of delta = 0)#
  # identification is repeated (shared peptide between two or more proteins)#

  uniqueHistoneForms <- NotIdentical |>
    dplyr::distinct(.data = _, across(!!select_exp) , .keep_all = TRUE)

  rlang::eval_tidy(uniqueHistoneForms)



# Iterate over the rows of uniqueHistoneForms
  for (i in seq_len(nrow(uniqueHistoneForms))) {


    PTM_value <- uniqueHistoneForms[[PTM_col]][i]
    var_value <- uniqueHistoneForms[[var_col]][i]

    if (!is.na(PTM_value)) { #in case non-modified peptide

      if(any(Duplicates[[PTM_col]] == PTM_value
             &
        Duplicates[[var_col]] == var_value)) {

        uniqueHistoneForms[[PTM_col]][i] <- paste0(PTM_value, "*")
     }
    }
  }



#count number of marked (*) PTMs
  num_marked_ptms = uniqueHistoneForms |>
    dplyr::filter(stringr::str_detect(!!PTM_col, "\\*")) |> nrow()



duration = round(as.numeric(Sys.time() - start_time), 1)
  cli::cli_alert_success('{num_marked_ptms} PTMs were marked (*) in {duration} s.')


  return(list(unique = uniqueHistoneForms, duplicate = Duplicates, freq_table = freq_table ))
}





