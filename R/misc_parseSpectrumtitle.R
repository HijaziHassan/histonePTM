misc_parseSpectrumtitle <- function(df, preprocess = c("MGFBoost", "MascotDistiller"), sheet = sheetName){

  preprocess = match.arg(preprocess)
  sheet = match.arg(sheet)


  column_of_interest <- "spectrum_title"

  if(!column_of_interest %in% colnames(df)) {

    cli::cli_alert_info("spectrum title column has either been already parsed or doesn't exist!")
  }

  else if (preprocess == "MGFBoost"){
    df <- df %>% separate(       #parse spectrum_title column into 3 sparate columns
      spectrum_title,
      into = c("cycle1", "cycle2", "scan", "scan_num", "rt_2", "rt_3", "file", "empty"),
      sep = ";", extra = "drop"
    ) %>%
      select(-c(cycle1, cycle2, scan, rt_2, rt_3, empty)) %>%
      mutate(file = as.character(str_split_i(file, ":", -1)),
             scan_num = parse_number(scan_num))



  }
  else if (preprocess == "MascotDistiller"){

    df <- df %>% separate(       #parse spectrum_title column into 3 sparate columns
      spectrum_title,
      into = c("n", "scan", "scan_num", "rt_2", "file"),
      sep = " "
    ) %>%
      select(-c(n, scan,  rt_2))%>%
      mutate(file = as.character(str_extract_all(file, "(?=\\w+\\d?+_\\d+).+(?=\\.)")))


  }

}
