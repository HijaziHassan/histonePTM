#' Histone Proteomics Data Analysis
#'
#' @description
#' Analyze histone bottom-up proteomics DDA data validated and quantified by `Proline` software.
#'
#' @param analysisfile `Proline` output excel file
#' @param metafile An excel file containing user-defined `SampleName` and `file` columns ('Condition', 'Bioreplicate' or 'TechReplicate' are `optional``).
#' @param hist_prot One of 4 histone proteins ('H3', 'H4', 'H2A', 'H2B'). If you want to analyze them all, choose "All".
#' @param labeling the labeling reagent used like 'PA' (default), 'TMA', 'PA_PIC', or 'none'. It's used to clear labeling modifications from PTM string by `misc_clearLabeling()`.
#' @param NA_threshold A cutt-off value of missing intensity values in which an identification is to be discarded (i.e. if 2, an ID with 2 or more NA will be discared).
#' @param extra_filter Either of 'no_me1', "K37un",  "no_me1_K37un", or "none" (default).
#' @param norm_method Normalization method. Either by 'peptide family' (default) or by "peptide_total". The latter depends on what is in your dataset and if you have prefiltered it or not.
#' 'no_me1' removes ALL peptides with unlabelled me1 . "no_me1_K37un" does the same but also removes H3K27-R40 peptides which are modified at K37.'none' does not do any filtration.
#' @param output_result Either `signle` or `multiple`. This will decided if all ids from different proteins are in one file (`single`) or in a separate file (`multiple`).
#'
#' @importFrom dplyr mutate filter across left_join any_of select starts_with arrange if_else rename where desc
#' @importFrom stringr str_detect str_split_i str_replace_all str_count str_trim
#' @importFrom tidyr nest drop_na
#' @importFrom writexl write_xlsx
#' @importFrom openxlsx writeData read.xlsx createWorkbook getSheetNames saveWorkbook addWorksheet
#' @importFrom rlang check_installed is_missing set_names
#' @importFrom purrr map map2 walk2 pluck
#' @importFrom cli cli_alert_warning cli_alert_success cli_abort cli_h1 cli_h2 cli_h3
#'
#' @return At least 3 excel files fragmented based on different filters:
#' \describe{
#'   \item{File1: fPSMs_analysisfile.xlsx}{
#'     \itemize{
#'       \item SHEET1: 'meta_data' — SampleName, file, and data columns (and Condition, Bioreplicate, TechReplicate if available).
#'       \item SHEET2: 'RawData' — The raw data as is with a subset of columns, some of which are renamed, and new columns like rt_diff, NA_count, etc., in addition to those parsed from the _spectrum_title_ column.
#'       \item SHEET3: 'Nt_IDs' — All peptides that are successfully N-terminally labelled. iRT are always included.
#'       \item SHEET4: '#NA<_NA_threshold' — IDs quantified with missing values less than the specified threshold.
#'       \item SHEET5: 'FullyMod_IDs' — All peptides that are N-terminally labelled and all Ks are modified endogenously or chemically.
#'       \item SHEET6: 'isob_coel_pep' — Peptides spotted by the `ptm_flagDupes()` function. Helps see which PTM × sequence combinations have a zero delta score.
#'       \item SHEET7: 'unique_IDs' — Unique IDs considered for quantification. Abundances are not normalized.
#'       \item SHEET8: 'unique_IDs_norm' — Same as the previous sheet but with normalized abundances.
#'       \item SHEET9: 'unique_IDs_norm_xxx' — Same as the previous sheet but normalized while discarding any peptide with what is defined by `extra_filter`.
#'     }
#'   }
#'   \item{File2: PTMsep_analysisfile.xlsx}{
#'     A file containing a summary of the identified PTMs. Each sheet per identified PTM is created.
#'     This file could be further fragmented into as many proteins as are present.
#'     To do this, change the argument `output_result` to 'multiple'.
#'   }
#'   \item{File3: hist_prot_analysisfile.xlsx}{
#'     A file per protein specified in `hist_prot`. Each file contains the identified peptide families, separated into sheets.
#'   }
#' }
#'
#'
#'
#'
#' @export
#'

analyzeHistone <- function(analysisfile,
                metafile,
                hist_prot= c('All','H3', 'H4', 'H2A', 'H2B'),
                labeling = c('PA', 'TMA', 'PIC_PA', "none"),
                NA_threshold,
                norm_method = c('peptide_family', 'peptide_total'),
                output_result= c('single', 'multiple'),
                extra_filter = c( "none", 'no_me1', "no_me1_K37un")){


  output_result = match.arg(output_result)
  hist_prot = match.arg(hist_prot)
  labeling = match.arg(labeling)
  norm_method = match.arg(norm_method)
  extra_filter = match.arg(extra_filter)
# Data 1st Check -----------------------------------------------

## packages -------
  rlang::check_installed(c("dplyr", 'stringr', 'purrr', 'tidyr', 'openxlsx', 'cli'))

## Proline excel output -------
  if(!file.exists(analysisfile)) cli::cli_abort("{analysisfile} was not found. Have you misspelled it?!")

## Excel spread sheet  -------
  sheetName= 'Best PSM from protein sets'

  excel_sheets_expected <- openxlsx::getSheetNames(file = analysisfile)
  if(!sheetName %in% excel_sheets_expected) cli::cli_abort("{sheetName} sheet is not found in {FileName}")

## User-defined rawfilenames -------
  if(!file.exists(metafile)) cli::cli_abort("{metafile} is not found or maybe misspelled!")


# create folders ------------------
  exp_folder_name  = tools::file_path_sans_ext(paste0("analysis/", analysisfile))
#img_folder_name  = tools::file_path_sans_ext(paste0("figures/", analysisfile))


# Overview ------------------
cli::cli({
  cli::cli_h1("Histone Analysis")
  cli::cli_h2("Sample: {exp_folder_name}")
  cli::cli_h3("Histones:  {hist_prot}")
  cli::cli_h3("Labeling method:  {labeling}")
})



 cat("\nCreating folders to store the results ...\n")
  misc_createFolder(foldername = "analysis")
  #misc_createFolder(foldername = "figures")
  misc_createFolder(foldername = exp_folder_name)
  #misc_createFolder(foldername = img_folder_name)



# Data 2nd Check -----------------------------------------------

## read files ---------
#readxl function seems to remove spaces upon reading, but to avoid loading extra packges, I used openxlsx
rawfilenames <- openxlsx::read.xlsx(xlsxFile = metafile) |>
    dplyr::mutate(dplyr::across(.cols = dplyr::where(is.character), .fns = ~stringr::str_trim(.x)))

proline_output <- openxlsx::read.xlsx(xlsxFile = analysisfile, sheet = sheetName) #read excel file into R


## extract metadata from Proline output ---------------

meta_names <- misc_extractMetaData(analysisfile)

## check matching

misc_checkFileNames(meta_names, rawfilenames, common_col = "file")

## create a named vector to store column names-------------
ColNames <- misc_useCustomNames(meta_names,rawfilenames) # this specific to Proline


## count number of samples ------------
numberofsamples <- sum(grepl("^abundance_", colnames(proline_output)))






# Code Execution ------------------------------------------------------------


## Data Prepare & Clean -----------------------------------------------

#### Process spectrum_title column------------------------

proline_output <- misc_parseSpectrumtitle(proline_output, col = "spectrum_title",
                                          scan_prefix = "first_scan:|Scan " , #to account for MGFBoost processing or Mascot Distiller
                                          file_prefix = "raw_file_identifier:|.*\\/" , sep = ";|\\.raw\\]") #here as well



#### merge real names with meta names and remove channel column ----------------
meta_names_merge <- dplyr::left_join(meta_names, rawfilenames, by = "file") |>
  dplyr::select(SampleName, file, dat, dplyr::any_of(c('Condition', 'BioReplicate', 'TechReplicate'))) |>
  tidyr::drop_na(SampleName, file, dat)

##### SHEET 1 -----------------------------------------------

SHEET1 = meta_names_merge

#add sample real names to the dataset
proline_output <- dplyr::left_join(proline_output, meta_names_merge,
                                   by = c("id_filename"="file"))



Histone <-  seq_getHistPeptide(df = proline_output, seq_col = sequence, histoneProtein = hist_prot)
#cli::cli_inform("{histoneProtein} histone peptides are extracted.")


#extract peptides and assign names to them (labelbyk is used here)

if ("All" %in% hist_prot) {
  prot_to_extract <- c('H3', 'H4', "H2A", 'H2B')
} else {
  prot_to_extract <- hist_prot
}


named <- purrr::map(.x = prot_to_extract,
             function(x) {rlang::set_names(purrr::pluck(sequenceDB, x, "labelbyK"),
                                           purrr::pluck(sequenceDB, x, "argC_frags")
             )
                                           }) |>  unlist()



add_irt <- rlang::set_names(purrr::pluck(sequenceDB, "iRT", "label"),
                   purrr::pluck(sequenceDB, "iRT", "tryp_frags")
)


named <- c(named, add_irt)


#give label to each sequence


Histone <- Histone |>
                #add column which labels each sequence by variant and K position
  dplyr::mutate(seq_stretch= stringr::str_replace_all(sequence, named),
                #extract Protein name without variants (H3 for H3 and H3.3, H2A for all H2A variants)
                protein = substr(stringr::str_split_i(seq_stretch, "_", 1), 1, 3) |>
                  #remove "." from H3.3, "m" from H3mm7/13 and "t" from "H3t"
                  gsub("[[:punct:]]|m|t", "", x= _)
  )



#select only columns of interest

Histone <- Histone |>
  dplyr::select(
  protein,
  sequence,
  seq_stretch,
  modifications,
  PTM = ptm_protein_positions,
  psm_score,
  ptm_score,
  moz = experimental_moz,
  rt_apex = master_elution_time,
  rt,
  z= charge,
  id_query = initial_query_id,
  id_scan,
  dat,
  id_filename,
  SampleName,
  dplyr::any_of(c(
    "Condition",
    "BioReplicate",
    'TechReplicate')),
  dplyr::starts_with(c("psm_count", "abundance_"))) |>
  dplyr::mutate(diff_rt = rt_apex - rt, .after =  rt)


#rename_PTMs
Histone <- Histone |>
  dplyr::mutate(
    PTM_proforma = ptm_toProForma(seq = sequence, mod= modifications, lookup = histptm_mass),
    PTM = ptm_beautify(PTM, software = 'Proline', lookup= histptm_lookup), .after = PTM) |>
  dplyr::select(-modifications) #no need for this column anymore.


#Count NA & PSMs rowwise

Histone <- Histone |>
  dplyr::mutate(
    NA_count = rowSums(is.na(dplyr::across(dplyr::starts_with("abundance_")))),
    PSM_count = rowSums(dplyr::across(dplyr::starts_with("psm_count")))
  )


#####SHEET 2 -----------------------------------------------
#raw data with only peptides based on user defined protein: H2A, H2B, H3, H4 or 'All". iRT is considered in all cases.
#renamed columns, renamed PTMs,
#new columns: 1) NA_count 2) PSM_count 3)PTM_proforma.
            # Parsed spectrum_title column to  4)'id_scan" and 5)"id_filename" columns.
            #6) "seq_stretch" which is label by K of each peptide.


SHEET2 <- Histone


#####SHEET 3 -----------------------------------------------
#only N-terminally labeled peptides (including miscleavages [if any])##
Nt_Histone <- Histone |>
  dplyr::filter(stringr::str_detect(PTM, "Nt") | stringr::str_detect(seq_stretch, "iRT"))

# if(any("iRT" %in% unique(Histone$variant))){
#   cli::cli_inform("Nt-terminally labeled peptides and iRTs are isolated.\n")}else{
#     cli::cli_inform("Nt-terminally labeled peptides are isolated.\n")}


SHEET3 <- Nt_Histone

#####SHEET 4 -----------------------------------------------
#Contain only identifications where number of NAs is below user-defined threshold
# if not defined, the threshold will be number of samples (i.e. only remove an id if it is not quantified.)

if(rlang::is_missing(NA_threshold)){NA_threshold = numberofsamples}

accepted_Histone <- Nt_Histone |>
  dplyr::filter(NA_count <= NA_threshold) |>
  dplyr::mutate(fully_modified = dplyr::if_else(stringr::str_detect(sequence, "K"),  #mark sequences having lysines in their sequences
                                                dplyr::if_else(
                                                  #for sequence containing K, mark those whose number of PTMs match number of Ks
                                                  stringr::str_count(sequence, "K") == stringr::str_count(PTM, "K"),
                                                  "TRUE",
                                                  "FALSE"),
                                                "Lysine_free_seq"))|>
  dplyr::mutate(fully_modified = dplyr::if_else(stringr::str_detect(seq_stretch, "iRT"), "iRT", fully_modified))

cli::cli_alert_warning("IDs with {NA_threshold} or more missing values were discarded.")

SHEET4 <- accepted_Histone


#####SHEET 5 -----------------------------------------------
#select peptides that are:  1) fully modified.i.e. all lysine bear modification (labeled or endo)
                          # 2)  are N-terminally modified and don't contain K.
                          # 3) iRTs
#Add new column: PTM_stripped: remove residues and their positions.


CompleteHistoneCases <- accepted_Histone |>
  #select iRTs, lysine-free sequences, and completely modified sequences.
  dplyr::filter(fully_modified %in% c("TRUE", "Lysine_free_seq", "iRT")) |>
  dplyr::mutate(
  PTM_stripped = ptm_beautify(PTM, software = 'Proline', lookup= histptm_lookup, residue = 'remove'),
  .after = PTM
  ) |>
  dplyr::arrange(dplyr::desc(psm_score))


SHEET5 <- CompleteHistoneCases


#####SHEET 6 & 7 -----------------------------------------------

#remove duplications where:
# the intensity value is identical for more than one ID (case of delta = 0)#
# identification is repeated (shared peptide between two or more proteins)#



dupesAnalysis <- ptm_flagDupes(df= CompleteHistoneCases,
                                    PTM_col = PTM,
                                    var_col = sequence,
                                    select_cols = dplyr::starts_with("abundance_"))


SHEET6 <- dupesAnalysis[[2]] # duplicate dataframe

SHEET7 <- uniqueHistoneForms <- dupesAnalysis[[1]] #unique-ID dataframe


#####SHEET 8 -----------------------------------------------

#same data as SHEET 7  but:
          #1) values are normalized
          #2) intensity columns are renamed.
          #4) Coef_Var is calculated per condition.
          #4) new column: PTM_unlabeled

SHEET8 <- uniqueHistoneForms |>
  quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                     grouping_var = sequence,
                     norm_method = norm_method) |>
  dplyr::rename(!!!dplyr::any_of(ColNames)) |>
  quant_coefVariation(df= _, df_meta= meta_names_merge,
                      int_col= dplyr::any_of(meta_names_merge$SampleName),
                      seq_col = sequence,
                      ptm_col = PTM,
                      format = 'wide') |>
  dplyr::mutate(PTM_unlabeled = misc_clearLabeling(PTM_stripped,
                                                   labeling = labeling),
 .after= 'PTM')


##names of the sheets where the previous data produced will be saved##
all_sheets <-
  c("meta_data",
    "RawData",
    "Nt_IDs",
    paste0("#NA<", NA_threshold),
    "FullyMod_IDs",
    "isob_coel_pep",
    "unique_IDs",
    "unique_IDs_norm"
  )

#####SHEET 9 -----------------------------------------------

if(extra_filter == "no_me1"){

  SHEET9 <- uniqueHistoneForms_nome1 <- uniqueHistoneForms |>
    dplyr::filter(!stringr::str_detect(PTM_stripped, "(?:[a-zA-Z]\\d*me1)\\b"))|>
    dplyr::mutate(PTM_unlabeled = misc_clearLabeling(PTM_stripped, labeling = labeling), .after= 'PTM') |>
    quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                       grouping_var = sequence) |>
    dplyr::rename(!!!dplyr::any_of(ColNames)) |>
    quant_coefVariation(df= _, df_meta= meta_names_merge,
                        int_col= dplyr::any_of(meta_names_merge$SampleName),
                        seq_col = sequence,
                        ptm_col = PTM,
                        format = 'wide')

  all_sheets <- c(all_sheets, "unique_IDs_norm_no_me1")
}else if(extra_filter == "K37un"){

  SHEET9 <-  uniqueHistoneForms |>
    dplyr::mutate(PTM_unlabeled = misc_clearLabeling(PTM_stripped, labeling = labeling), .after= 'PTM') |>
    dplyr::filter(
      !stringr::str_detect(sequence, 'KS.P..GGVKKPHR') |
        (stringr::str_detect(sequence, 'KS.P..GGVKKPHR') & stringr::str_detect(PTM_unlabeled, '^.+-.+-un$'))
    ) |>
    quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                       grouping_var = sequence) |>
    dplyr::rename(!!!dplyr::any_of(ColNames)) |>
    quant_coefVariation(df= _, df_meta= meta_names_merge,
                        int_col= dplyr::any_of(meta_names_merge$SampleName),
                        seq_col = sequence,
                        ptm_col = PTM,
                        format = 'wide')
  all_sheets <- c(all_sheets, "unique_IDs_norm_K37un")

}else if(extra_filter == "no_me1_K37un"){

  SHEET9 <- uniqueHistoneForms_nome1 <- uniqueHistoneForms |>
    dplyr::mutate(PTM_unlabeled = misc_clearLabeling(PTM_stripped, labeling = labeling), .after= 'PTM') |>
    dplyr::filter(!stringr::str_detect(PTM_stripped, "(?:[a-zA-Z]\\d*me1)\\b"),
                  !stringr::str_detect(sequence, 'KS.P..GGVKKPHR') |
                    (stringr::str_detect(sequence, 'KS.P..GGVKKPHR') & stringr::str_detect(PTM_unlabeled, '^.+-.+-un$')))|>
    quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                       grouping_var = sequence) |>
    dplyr::rename(!!!dplyr::any_of(ColNames)) |>
    quant_coefVariation(df= _, df_meta= meta_names_merge,
                        int_col= dplyr::any_of(meta_names_merge$SampleName),
                        seq_col = sequence,
                        ptm_col = PTM,
                        format = 'wide')
all_sheets <- c(all_sheets, "unique_IDs_norm_no_me1_K37un")
}else if(extra_filter == none){

  all_sheets <- all_sheets[-length(all_sheets)]
  SHEET9 = NULL

}

list_dfs= list(SHEET1, SHEET2, SHEET3, SHEET4, SHEET5, SHEET6, SHEET7, SHEET8, SHEET9)
list_dfs <- Filter(Negate(is.null), list_dfs)

#setNames(list_dfs, all_sheets)


## Output Files -----------------------------------


##### --------- output file 1: Proteins ---------#######


if(output_result == "multiple"){
#filtered_dfs <- list()

for (prot in prot_to_extract) {

  list_prot <- list()

  for (i in seq_along(list_dfs)) {

    df <- list_dfs[[i]]

    if ('protein' %in% colnames(df)) { # to skip metadata

      # df_filtered <- df |>
      #   dplyr::filter(protein == prot)

      list_prot[[i]] <- df[df$protein == prot,]

    } else {

      list_prot[[i]] <- df
    }
  }



  ##create a workbook to save all the datatframes that will be produced##
  wb = openxlsx::createWorkbook()

  ##Add the sheets to a workbook (excel file)##
  purrr::map(.x = all_sheets , .f = ~openxlsx::addWorksheet(wb=wb, sheetName = .x))

  # Add each filtered data frame to a separate sheet
  purrr::map2(
    .x = list_prot,
    .y = all_sheets,
    .f = ~ openxlsx::writeData(wb = wb, x = .x, sheet = .y)
  )

  ##save the file containing the ##cleaned data ##
  folderName= tools::file_path_sans_ext(analysisfile)
  openxlsx::saveWorkbook(wb=wb,
                         file = paste0("./analysis/",folderName,"/",  paste0(prot, "_fPSMs_", analysisfile)),
                         overwrite = TRUE)

  cli::cli_alert_success('An excel file for {prot} histone protein is created.')


}}else if(output_result == "single"){


##create a workbook to save all the datatframes that will be produced##
wb = openxlsx::createWorkbook()

##Add the sheets to a workbook (excel file)##
purrr::map(.x = all_sheets , .f = ~openxlsx::addWorksheet(wb=wb, sheetName = .x))

# Add each filtered data frame to a separate sheet
purrr::map2(
  .x = list_dfs,
  .y = all_sheets,
  .f = ~ openxlsx::writeData(wb = wb, x = .x, sheet = .y)
)


##save the file containing the ##cleaned data ##
folderName= tools::file_path_sans_ext(analysisfile)
openxlsx::saveWorkbook(wb=wb,
                       file = paste0("./analysis/",folderName,"/fPSMs_", analysisfile),
                       overwrite = TRUE)

cli::cli_alert_success('An excel file containing cleaned and renamed data of {prot_to_extract} is created.')

}


##### --------- output file 2: peptides ---------#######

#determine the sheet which will be used for the ptm-centered and peptide-centered excel files to be generated.

df_tobe_splitted <- switch(extra_filter,
             "none" = SHEET8,
             "no_me1" = SHEET9,
             "K37un"= SHEET9,
             "no_me1_K37un" = SHEET9,
             SHEET8
)

df_tobe_splitted |>

    tidyr::nest(.by = protein) |>
  dplyr::mutate(
    split_data = purrr::map(data, ~ split(.x, .x['seq_stretch'])),
    # named dfs will be saved as sheets with write_xlsx
    save_files = purrr::walk2(.x = protein,
                              .y = split_data,
                             ~ writexl::write_xlsx(.y, paste0('./analysis/', folderName,"/", .x, ".xlsx")))
  )


for(protein in prot_to_extract){
cli::cli_alert_success('{protein}.xlsx is saved.')
}

##### --------- output file 3: PTMs---------#######
#create a sheet per PTM#



wb_ptm = openxlsx::createWorkbook()
ident_ptms <- strsplit(df_tobe_splitted$PTM_unlabeled, "-") |> unlist() |> unique()
purrr::map(.x = ident_ptms , .f = ~openxlsx::addWorksheet(wb=wb_ptm, sheetName = .x))

for(ptm in ident_ptms){

  df_tobe_splitted |>
    dplyr::filter(stringr::str_detect(PTM_unlabeled, ptm)) |>
    openxlsx::writeData(wb = wb_ptm, x = _, sheet = ptm)
}


openxlsx::saveWorkbook(wb=wb_ptm,
             file = paste0('./analysis/', folderName, "/PTMsep_", analysisfile),
             overwrite = TRUE)


cli::cli_alert_success('An excel file summarizing the IDs per each {ident_ptms} is created.')


cat("\n", date(), "\n")
cat(" A plus dans le bus ^_^!")


  }



