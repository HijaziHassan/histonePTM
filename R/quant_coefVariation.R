#' Calculate Coefficient of Variation (CV) for peptide intensities
#'
#' Computes mean, standard deviation, and coefficient of variation (CV)
#' for each peptide (and PTM) across replicate samples within each condition.
#' Supports both wide and long data formats.
#'
#' @param df A data frame containing peptide intensities.
#' @param df_meta Metadata with `SampleName` and `Condition` columns.
#'   Required for `format = "wide"`.
#' @param int_col Character vector of intensity column names, or a single
#'   intensity column name for long format.
#' @param seq_col Unquoted column name for sequence identifier.
#' @param ptm_col Unquoted column name for PTM information (optional).
#' @param format Either `"wide"` (default) or `"long"`.
#' @param scale Multiplier for CV. Default `100` (percent).
#' @param min_replicates Minimum replicates required. Default `2`.
#' @param ... Additional grouping columns for long format.
#'
#' @return The input data frame with added columns:
#'   `avg_<condition>`, `sd_<condition>`, `cv_<condition>`.
#'
#' @importFrom rlang is_missing sym
#' @importFrom dplyr select summarise filter group_by left_join all_of mutate
#' @importFrom cli cli_abort cli_alert_warning
#' @importFrom stats sd
#' @export
quant_coefVariation <- function(df,
                                df_meta = NULL,
                                int_col,
                                seq_col,
                                ptm_col = NULL,
                                format = c("wide", "long"),
                                scale = 100,
                                min_replicates = 2,
                                ...) {

  format <- match.arg(format)

  # ---- Input validation ----
  if (rlang::is_missing(df)) {
    cli::cli_abort("Did you forget to pass your dataset to the `df` argument?")
  }
  if (rlang::is_missing(seq_col)) {
    cli::cli_abort("You did not provide the name of the sequence column in `seq_col`.")
  }

  # ---- Empty data guard ----
  if (nrow(df) == 0) {
    cli::cli_alert_warning("Input data frame is empty. Returning as is.")
    return(df)
  }

  # ---- Wide format ----
  if (format == "wide") {
    if (rlang::is_missing(df_meta)) {
      cli::cli_abort(c(
        "If data is in `format = 'wide'`, `df_meta` must be provided.",
        "x" = "`df_meta` argument is missing.",
        "i" = "Provide a data frame with `SampleName` and `Condition` columns."
      ))
    }

    # Get intensity column names
    int_cols <- df |> dplyr::select(dplyr::all_of({{ int_col }})) |> names()

    if (rlang::is_empty(int_cols)) {
      cli::cli_abort(c(
        "Missing or incorrect input.",
        "x" = "The intensity columns' names were not found.",
        "i" = "Check the names provided in `int_col`."
      ))
    }

    # Clean metadata
    df_meta <- df_meta[, !grepl("^\\.{3}[1-9]$|^$", colnames(df_meta))]
    df_meta <- df_meta |> dplyr::filter(SampleName %in% int_cols)

    if (nrow(df_meta) == 0) {
      cli::cli_abort("No matching samples found between `df` and `df_meta`.")
    }

    # Map samples to conditions
    sample_cond <- split(df_meta$SampleName, df_meta$Condition)
    unique_cond <- unique(df_meta$Condition)
    unique_cond <- unique_cond[!is.na(unique_cond)]

    df_cv <- df

    for (cond in unique_cond) {
      indices <- which(colnames(df_cv) %in% intersect(int_cols, sample_cond[[cond]]))

      if (length(indices) < min_replicates) {
        cli::cli_alert_warning(
          "Condition '{cond}' has only {length(indices)} replicate(s). CV will be NA."
        )
      }

      avg_col <- paste0("avg_", cond)
      sd_col  <- paste0("sd_", cond)
      cv_col  <- paste0("cv_", cond)

      df_cv <- df_cv |>
        dplyr::mutate(
          !!avg_col := rowMeans(df_cv[, indices, drop = FALSE], na.rm = TRUE),
          !!sd_col  := apply(df_cv[, indices, drop = FALSE], 1, sd, na.rm = TRUE)
        )

      # Compute CV safely (avoid division by zero)
      df_cv <- df_cv |>
        dplyr::mutate(
          !!cv_col := dplyr::if_else(
            .data[[avg_col]] == 0 | is.na(.data[[avg_col]]) | length(indices) < min_replicates,
            NA_real_,
            round((.data[[sd_col]] / .data[[avg_col]]) * scale, 2)
          )
        )
    }

    return(df_cv)

    # ---- Long format ----
  } else if (format == "long") {

    # Capture additional grouping columns
    extra_groups <- rlang::enquos(...)

    if ("Condition" %in% colnames(df)) {
      df_cv <- df |>
        dplyr::select({{ seq_col }}, {{ ptm_col }}, dplyr::all_of(int_col), !!!extra_groups) |>
        dplyr::group_by(Condition, {{ ptm_col }}, {{ seq_col }}, !!!extra_groups) |>
        dplyr::summarise(
          n_replicates = dplyr::n(),
          mean_val = mean(.data[[int_col]], na.rm = TRUE),
          sd_val   = sd(.data[[int_col]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          CV = dplyr::if_else(
            n_replicates >= min_replicates & mean_val != 0 & !is.na(mean_val),
            round((sd_val / mean_val) * scale, 2),
            NA_real_
          )
        ) |>
        dplyr::arrange(CV)
    } else if (rlang::is_missing(df_meta)) {
      cli::cli_abort(c(
        "If data is in `format = 'long'` and contains no `Condition` column,",
        "`df_meta` must be provided.",
        "x" = "`Condition` column is missing.",
        "i" = "Provide `df_meta` with `SampleName` and `Condition` columns."
      ))
    } else {
      df_cv <- df |>
        dplyr::select({{ seq_col }}, {{ ptm_col }}, dplyr::all_of(int_col), !!!extra_groups) |>
        dplyr::left_join(df_meta, by = "SampleName") |>
        dplyr::group_by(Condition, {{ ptm_col }}, {{ seq_col }}, !!!extra_groups) |>
        dplyr::summarise(
          n_replicates = dplyr::n(),
          mean_val = mean(.data[[int_col]], na.rm = TRUE),
          sd_val   = sd(.data[[int_col]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          CV = dplyr::if_else(
            n_replicates >= min_replicates & mean_val != 0 & !is.na(mean_val),
            round((sd_val / mean_val) * scale, 2),
            NA_real_
          )
        ) |>
        dplyr::arrange(CV)
    }

    # Warn if any CVs are NA
    if (any(is.na(df_cv$CV))) {
      cli::cli_alert_warning(
        "Some CV values are NA (insufficient replicates or mean = 0)."
      )
    }

    return(df_cv)
  }
}
