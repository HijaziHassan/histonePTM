
#create costum names to replace dat file coded abundance and psm count headers
misc_useCostumNames <- function(raw_df, user_df){
  # the unmatced argument is not working.
  abn_ColNames <- dplyr::left_join(raw_df,
                                   user_df,
                                   by= "file",
                                   keep = FALSE,
                                   unmatched = "drop") |>
    mutate(new_name = paste0("abundance_", channel),
           .keep = "unused")


  psmcount_ColNames <- left_join(raw_df,
                                 user_df ,
                                 by= "file",
                                 keep = FALSE,
                                 unmatched = "drop") |>
    mutate(SampleName = paste0(SampleName, "_psmCount"),
           new_name = paste0("psm_count_", channel),
           .keep = "unused")

  ColNames <- bind_rows(abn_ColNames, psmcount_ColNames) %>%
    select(-c(file, dat, any_of(c('Condition', 'Bioreplicate')))) %>%
    deframe()

  ColNames
}
