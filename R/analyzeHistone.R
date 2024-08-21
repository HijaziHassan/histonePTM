#' Histone Proteomics Data Analysis
#'
#' @description
#' Analyze histone bottom-up proteomics DDA data validated and quantified by `Proline` software.
#'
#' @param analysisfile `Proline` output excel file
#' @param metafile An excel file containing user-defined `SampleName` and `file` columns ('Condition', 'Bioreplicate' or 'TechReplicate' are `optional``).
#' @param hist_prot One of 4 histone proteins ('H3', 'H4', 'H2A', 'H2B'). If you to analyze them all use "All".
#' @param NA_threshold A filter value below which an identification having this value of missing intensity value(s) or more is discarded.
#' @param output_result Either `signle` or `multiple`. This will decided if all ids from different proteins are in one file (`single`) or in a separate file (`multiple`).
#'
#' @import dplyr
#' @importFrom stringr str_detect str_split_i str_replace_all str_count
#' @importFrom tidyr nest drop_na
#' @importFrom openxlsx writeData read.xlsx createWorkbook getSheetNames saveWorkbook addWorksheet
#' @importFrom rlang check_installed is_missing set_names
#' @importFrom purrr map map2 walk2 pluck
#' @importFrom cli cli_alert_warning cli_alert_success cli_abort cli_h1 cli_h2 cli_h3
#' @return At least 3 excel files.
#' @export
#'

analyzeHistone <- function(analysisfile,
                metafile,
                hist_prot,
                NA_threshold,
                output_result= c('single', 'multiple')){


# Data 1st Check -----------------------------------------------

## packages -------
  rlang::check_installed(c("dplyr", 'stringr', 'purrr', 'tidyr', 'openxlsx', 'cli'))

## Proline excel output -------
  if(!file.exists(analysisfile)) cli::cli_abort("{analysisfile} is not found or maybe misspelled!")

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
})



 cat("\nCreating folders to store the results ...\n")
  misc_createFolder(foldername = "analysis")
  #misc_createFolder(foldername = "figures")
  misc_createFolder(foldername = exp_folder_name)
  #misc_createFolder(foldername = img_folder_name)



# Data 2nd Check -----------------------------------------------

## read files ---------

rawfilenames <- openxlsx::read.xlsx(xlsxFile = metafile)

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
                                          scan_prefix = "first_scan:" ,
                                          file_prefix = "raw_file_identifier:" , sep = ";")



#### merge real names with meta names and remove channel column ----------------
meta_names_merge <- dplyr::left_join(rawfilenames, meta_names, by = "file") |>
  dplyr::select(SampleName, file, dat, any_of(c('Condition', 'Bioreplicate'))) |>
  tidyr::drop_na()

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


named <- purrr:::map(.x = prot_to_extract,
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
                  #remove "." from H3.3 and "m" from H3mm7/13
                  gsub("[[:punct:]]|m", "", x= _)
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
    PTM_proforma = ptm_toProForma(seq = sequence, mod= modifications),
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
  dplyr::mutate(PTM_stripped1 = stringr::str_replace_all(
    PTM, # rename PTM for simiplicity (remove Nt- and Kxx)
    c(".+Nt-" = "",
      ".+Nt" = "Nt_noK",
      "K\\d+" = "")
  ),
  PTM_stripped2 = ptm_beautify(PTM, software = 'Proline', lookup= histptm_lookup, residue = 'remove'),
  .after = PTM
  ) |>
  dplyr::arrange(desc(psm_score))


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

SHEET8 <- uniqueHistoneForms |>
  quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                     grouping_var = sequence) |>
  dplyr::rename(!!!dplyr::any_of(ColNames)) |>
  quant_coefVariation(df= _, df_meta= rawfilenames,
                      int_col= dplyr::any_of(rawfilenames$SampleName),
                      seq_col = sequence,
                      ptm_col = PTM,
                      format = 'wide')

#####SHEET 9 -----------------------------------------------

uniqueHistoneForms_nome1 <- uniqueHistoneForms |>
  dplyr::filter(!stringr::str_detect(PTM_stripped2, "(?:[:upper:]*\\d*me1)\\b")) |>
  dplyr::mutate(PTM_unlabeled = misc_clearLabeling(PTM_stripped2), .after= 'PTM') |>
  quant_relIntensity(select_cols = dplyr::starts_with("abundance_"),
                     grouping_var = sequence) |>
  dplyr::rename(!!!dplyr::any_of(ColNames)) |>
  quant_coefVariation(df= _, df_meta= rawfilenames,
                      int_col= dplyr::any_of(rawfilenames$SampleName),
                      seq_col = sequence,
                      ptm_col = PTM,
                      format = 'wide')


SHEET9 <- uniqueHistoneForms_nome1



##names of the sheets where the previous data produced will be saved##
all_sheets <-
  c("meta_data",
    "RawData",
    "Nt_IDs",
    paste0("#NA<", NA_threshold),
    "FullyMod_IDs",
    "isob_coel_pep",
    "unique_IDs",
    "unique_IDs_norm",
    "unique_IDs_norm_nome1"
  )





list_dfs= list(SHEET1, SHEET2, SHEET3, SHEET4, SHEET5, SHEET6, SHEET7, SHEET8, SHEET9)
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

# print(prot_to_extract)
# seq_spreadIntoSheets()

#create a sheet per variant and save the separated data
#extract_variants(uniqueHistoneForms_normalized, histoneProtein =  histoneProtein)

suppressWarnings({SHEET9 |>

    tidyr::nest(.by = protein) |>
  dplyr::mutate(
    split_data = purrr::map(data, ~ split(.x, .x['seq_stretch'])),
    # named dfs will be saved as sheets with write_xlsx
    save_files = purrr::walk2(.x = protein,
                             .y = split_data,
                             ~ writexl::write_xlsx(.y, paste0('./analysis/', folderName,"/", .x, ".xlsx")))
  )

})

cli::cli_alert_success('An excel file per {prot_to_extract} is created separately.')

##### --------- output file 3: PTMs---------#######
#create a sheet per PTM#



wb_ptm = openxlsx::createWorkbook()
ident_ptms <- strsplit(SHEET9$PTM_unlabeled, "-") |> unlist() |> unique()
purrr::map(.x = ident_ptms , .f = ~openxlsx::addWorksheet(wb=wb_ptm, sheetName = .x))

for(ptm in ident_ptms){

  SHEET9 |>
    dplyr::filter(stringr::str_detect(PTM_unlabeled, ptm)) |>
    openxlsx::writeData(wb = wb_ptm, x = _, sheet = ptm)
}


openxlsx::saveWorkbook(wb=wb_ptm,
             file = paste0('./analysis/', folderName, "/PTMsep_", analysisfile),
             overwrite = TRUE)


cli::cli_alert_success('An excel file summarizing ids per each {ident_ptms} is created.')


cat("\n", date(), "\n")
cat(" A plus dans le bus ^_^!")


  }



