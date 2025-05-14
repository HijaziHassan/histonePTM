#' Correct of Naturally Occuring Abundances
#' @description
#' Subtract naturally occuring abundances of MS XICs in stable isotope tracing experiments
#'
#'
#' @param df A dataframe as PICor input.
#' @param time_col The time column.
#' @param t0 Value of the reference sample in the time column ( e.g. "0" , "t0", "T0" ...)
#' @param cond_col Condition column (Concentration, disease vs non-disease ...)
#' @param avg_by (optional) A column to average by (e.g.TechReplicate).
#'
#' @importFrom readr parse_number
#' @importFrom dplyr group_by group_split group_keys pull select bind_rows bind_cols rowwise mutate filter across matches all_of summarise if_else
#' @importFrom purrr map
#' @importFrom stringr str_to_lower str_glue str_flatten
#' @importFrom tibble deframe
#' @importFrom tidyr pivot_longer
#' @importFrom cli cli_abort
#' @importFrom rlang sym
#' @return A dataframe with corrected and normalized values.
#' @export
#'
#'
quant_manIsoCorreciton <- function(df, time_col, t0, cond_col, avg_by) {


  # Step 1: Extract columns of isotopes
  col_names <- colnames(df)
  ## Step 1.1 classify them into odd and even
  all_isotopes <- classify_isotopes(col_names)
  ## Step 1.2: isolate the even ones (acetyl has two heavy carbons)
  even_iso <- all_isotopes[['even_iso']]

  ## Step 1.3: get-rid of odd C13 istopes
  df <- df |>
    dplyr::select(!dplyr::matches(all_isotopes[['odd_iso']]))

  ### Step 1.3.1
  if(!missing(avg_by)){

    df <- get_avg(df, col = {{avg_by}})
  }

  #Step 2: Split the dataframe by `cond_col` into a named list

  split_dfs <- df |>
    dplyr::group_by({{cond_col}}) |>
    dplyr::group_split() |>
    setNames(df |>
               dplyr::group_by({{cond_col}}) |>
               dplyr::group_keys() |>
               dplyr::pull({{cond_col}})
    )


  # Step 3: Correct & Normalize
  split_dfs <- purrr::map(names(split_dfs), function(name) {

    ## Step 3.1: grab one df at each go
    df_splitted <- split_dfs[[name]]

    ## Step 3.2: Calculate M+0 correction ratios for t0
    ratio <- M0_corr_ratio(df= df_splitted, time_col= {{time_col}}, t0 = t0, even_iso = even_iso)

    ## Step 3.3 : Subtract naturally occuring abundances
    df_corrected <- correct_isotopes(df= df_splitted , time_col= {{time_col}}, t0 = t0, M0_ratios = ratio, even_iso= even_iso)

    ## Step 3.4: Normalize corrected isotopic distributions
    df_norm <- perc_incorporation(df_corrected)

    return(df_norm)

  }) |> dplyr::bind_rows()

  split_dfs
}




# Supporting Functions -----------------------------------------

#' @noRd
classify_isotopes <- function(col_names) {
  #This is based on the output of PiCorr software
  c13_isotopes <- grep("C13", col_names, value = TRUE)
  if (length(c13_isotopes) == 0) {
    warning("No C13 isotopes detected in the provided column names.")
    return(list(even_iso = character(), odd_iso = character()))
  }

  # Parse numeric parts of the isotope labels
  isotope_numbers <- readr::parse_number(c13_isotopes)

  even_isotopes <- c13_isotopes[!is.na(isotope_numbers) & isotope_numbers %% 2 == 0]
  odd_isotopes <- c13_isotopes[!is.na(isotope_numbers) & isotope_numbers %% 2 != 0]


  return(list(even_iso = even_isotopes, odd_iso = odd_isotopes))
}

#' @noRd
M0_corr_ratio <- function(df, time_col, t0, even_iso) {

  if (length(even_iso) == 0L) {
    cli::cli_abort(c(
      'x' = 'No columns of isotopes intensities were detected.',
      'i' = 'Are they represented as follows: "2C13, 4C13, ..."?')
    )
  }

  purrr::map(even_iso, ~{
    df |>
      dplyr::filter(stringr::str_to_lower({{time_col}}) == stringr::str_to_lower(t0)) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        "ratio_{.x}" := !!rlang::sym(.x) / `No label`,
        .keep = "none"
      ) |>
      dplyr::ungroup()
  }) |>
    dplyr::bind_cols() |>
    tidyr::pivot_longer(everything(),
                        names_to = 'isotope',
                        values_to = 'ratio') |>
    tibble::deframe()
}

#' @noRd
correct_isotopes <- function(df, time_col, t0, M0_ratios, even_iso) {

  if (length(even_iso) == 0L) {
    cli::cli_abort("No even isotopes were provided. Ensure `even_iso` is a valid character vector.")
  }

  if (length(M0_ratios) == 0L) {
    cli::cli_abort("The ratio vector is empty. Ensure `ratio` is a valid named vector.")
  }

  # M+0 is only C12
  df_corr <- df |>
    dplyr::mutate(`M+0` = `No label`)

  # Parse numeric parts of even isotopes
  iso_dist <- readr::parse_number(even_iso)

  # Correct for each even isotope dynamically
  for (i in seq_along(iso_dist)) {
    iso_num <- iso_dist[i]
    col <- even_iso[i]
    ratio_name <- paste0("ratio_", col)

    # Check if ratio exists
    if (!ratio_name %in% names(M0_ratios)) {
      cli::cli_abort("Ratio '{ratio_name}' not found in {.arg M0_ratios}.")
    }

    # Construct the correction formula dynamically
    if (i == 1) {
      # First isotope (M+2): Only subtract `No label`
      correction_formula <- stringr::str_glue("`{col}` - `No label` * M0_ratios[['{ratio_name}']]")
    } else {
      # Subsequent isotopes (M+4, M+6, ...): Include contributions from previous M+
      prev_iso_dist <- iso_dist[1:(i - 1)]
      correction_terms <- stringr::str_flatten(
        stringr::str_glue("`M+{prev_iso_dist}` * M0_ratios[['ratio_{rev(prev_iso_dist)}C13']]"),
        collapse = " - "
      )
      correction_formula <- stringr::str_glue("`{col}` - `No label` * M0_ratios[['{ratio_name}']] - {correction_terms}")
    }

    # Apply the correction dynamically
    df_corr <- df_corr |>
      dplyr::mutate(
        !!rlang::sym(paste0("M+", iso_num)) := dplyr::if_else(
          stringr::str_to_lower({{time_col}}) == t0,
          0,
          eval(parse(text = correction_formula))
        )
      )
  }

  return(df_corr)
}

#' @noRd
# normalize by sum of the whole (even) isotopic profile
perc_incorporation <- function(df_corrected) {

  df_normalized <- df_corrected |>
    dplyr::mutate(dplyr::across(dplyr::matches('M\\+\\d+'),
                                ~ ( ./ rowSums(dplyr::across(dplyr::matches('M\\+\\d+')), na.rm = TRUE)
                                ) * 100))

  return(df_normalized)
}


#' @noRd
#Avg intensities
get_avg <- function(df, col){
  df |>
    dplyr::summarise(dplyr::across(c(dplyr::all_of(c("No label")),
                                     dplyr::all_of(c(dplyr::ends_with("C13")))),
                                   ~ mean(.x, na.rm = TRUE)),
                     .by = {{col}})

}
