#' Separate peptides into spreadsheets
#'
#' @description
#' Each sequence/variant will be saved to a separate excel worksheet. This facilitate focus on a specific site or peptide.
#'
#' @param df A dataframe
#' @param histoneProtein Three options 'H3-H4', 'H3' or 'H4'. In none is selected 'H3-H4' will be selected by default.
#' @param split_by name of the column containing sequence labels. The values will be the name of the excel sheets.
#' @param label_type label type to match by from \code{SequenceDB}. Two options 'labelbyK' (default) and 'label'.
#' @param FileName In case you want to save the file as \code{.xlsx}, provide a name with this extension.
#'
#' @return A list of \code{Workbook} object(s), each with multiple excel worksheets.
#' @export
#' @import openxlsx
#' @importFrom tidytable group_split
#' @importFrom dplyr filter
#' @importFrom purrr walk walk2 pluck
#'
seq_spreadIntoSheets <- function(df,
                                 split_by= "variant",
                                 histoneProtein = c('H3-H4', 'H3', 'H4'),
                                 label_type = c('labelbyK', 'label'),
                                 FileName =""){

  histoneProtein = match.arg(histoneProtein)
  split_by = match.arg(split_by)
  label_type = match.arg(label_type)




  if(histoneProtein == "H3-H4"){

    wb_varaintsH3 = openxlsx::createWorkbook()
    wb_varaintsH4 = openxlsx::createWorkbook()

    dfH3 <- df |>
      dplyr::filter(sequence %in% purrr::pluck(sequenceDB, "H3", "argC_frags"))


    dfH4 <- df |>
      dplyr::filter(sequence %in% purrr::pluck(sequenceDB, "H4", "argC_frags"))

    #split data according to variants (H3canon, H3.3, H3mm7, and H3mm13, ....)
    variant_splitsH3 <- tidytable::group_split(dfH3, split_by, .named = TRUE)
    variant_splitsH3 <- variant_splitsH3[order(match(names(variant_splitsH3), purrr::pluck(sequenceDB, "H3", "labelbyK")))]

    variant_splitsH4 <- tidytable::group_split(dfH4, split_by, .named = TRUE)
    variant_splitsH4 <- variant_splitsH4[order(match(names(variant_splitsH4), purrr::pluck(sequenceDB, "H4", "labelbyK")))]

    #Extract variants' names to set them as sheet titles.
    split_namesH3 <- unique(dfH3[[split_by]])
    split_namesH3 <- split_namesH3[order(match(split_namesH3, purrr::pluck(sequenceDB, "H3", label_type)))]

    split_namesH4 <- unique(dfH4[[split_by]])
    split_namesH4 <- split_namesH4[order(match(split_namesH4, purrr::pluck(sequenceDB, "H4", label_type)))]


    ##store each split group/varait in its corresponding sheet matching its name##

    #H3
    split_namesH3 |>
      purrr::walk( ~ openxlsx::addWorksheet(wb_varaintsH3, sheetName = .x))

    purrr::walk2(.x = split_namesH3,
                 .y = variant_splitsH3,
                 .f = ~ openxlsx::writeData(wb = wb_varaintsH3, sheet = .x, x = .y),
                 .progress = TRUE)
    #H4
    split_namesH4 |>
      purrr::walk( ~ openxlsx::addWorksheet(wb_varaintsH4, sheetName = .x))

    purrr::walk2(.x= split_namesH4,
          .y= variant_splitsH4,
          .f = ~ openxlsx::writeData(wb = wb_varaintsH4, sheet = .x, x = .y),
          .progress = TRUE)



    if(FileName !=""){

        ##save the file containing ##data per variant##

    fileH3 = paste0('H3seq_', FileName)
    fileH4 = paste0('H4seq_', FileName)

    openxlsx::saveWorkbook(wb=wb_varaintsH3, file = fileH3 , overwrite = TRUE)
    openxlsx::saveWorkbook(wb=wb_varaintsH4, file = fileH4, overwrite = TRUE)

    cat("varaints/sequences are separated.\nFile", "\"", basename(fileH3), "\"", "is saved.\nFile", "\"", basename(fileH4), "\"", "is saved.\n\n")

    }
    return(list(wb_varaintsH3, wb_varaintsH4))

  }

  else if (histoneProtein == "H3" | histoneProtein == "H4"){

    # vector containing all H3 and H4 peptides
    peptide_order <- c(purrr::pluck(sequenceDB, "H3", label_type), purrr::pluck(sequenceDB, "H4", label_type))
    wb_varaints = openxlsx::createWorkbook()

    #split data according to variants (H3canon, H3.3, H3mm7, and H3mm13, ....)
    variant_splits <- tidytable::group_split(df, split_by, .named = TRUE)
    variant_splits <- variant_splits[order(match(names(variant_splits), peptide_order))]

    #Extract variants' names to set them as sheet titles.
    split_names <- unique(df[[split_by]])
    split_names <- split_names[order(match(split_names, peptide_order))]



    ##store each split group/varait in its corresponding sheet matching its name##
    split_names |>
      purrr::walk( ~ openxlsx::addWorksheet(wb_varaints, sheetName = .x))

    purrr::walk2(.x= split_names,
                 .y= variant_splits,
                 .f = ~ openxlsx::writeData(wb = wb_varaints, sheet = .x, x = .y), .progress = TRUE)



    if(FileName !=""){

    filePSM = paste0('Histseq_', histoneProtein,'_', FileName)
    # ##save the file containing ##data per variant##
    openxlsx::saveWorkbook(wb = wb_varaints, file = filePSM, overwrite = TRUE)

    cli::cli_alert_success("varaints/sequences are separated. File '{filePSM}' is saved.\n\n")
    }

    return(list(wb_varaints))

  }

}
