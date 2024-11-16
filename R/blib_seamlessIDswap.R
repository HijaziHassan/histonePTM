#main function
#' Change best spectrum ID in \emph{.blib} spectral library
#'
#' @description
#' To swap IDs between redundant library and main library in case of wrong peak picking or (near)-coeluting isobaric peptides
#' in Skyline software.
#'
#'
#' @param db_main The name of the main Skyline \code{.blib} library file (\code{string}).
#' @param db_redundant The name of the redundant Skyline \code{.blib} library file (\code{string}).
#' @param rt Retention time (\code{numeric}) of the ID to be set as the best spectrum (rounded to first decimal).
#' Also accepts numeric vector.
#' @param mz \emph{m/z} (\code{numeric}) value of the target ID (it will be truncated to two decimal places).
#' Also accepts numeric vector.
#' @param tol The tolerance (\code{optional}) to be used to match the provided \emph{m/z} (\code{numeric}).
#' The default value is used to have a perfect match.
#' @param file MS raw file name without extension (\code{optional}). In case of multiple matches to \code{rt},
#'  it is necessary to provide the \code{file} name to force specific ID selection (the file name is found at the top of
#' chromatogram window of each sample in Skyline interface).
#' Also accepts character vector.
#'
#' @importFrom dplyr select filter mutate collect tbl case_when left_join join_by pull bind_rows summarise
#' @importFrom cli cli_abort cli_alert_warning cli_alert_danger
#' @importFrom purrr pmap
#' @importFrom DBI dbRollback dbConnect dbDisconnect dbExecute dbBegin dbWriteTable dbCommit
#' @importFrom RSQLite SQLite
#'
#' @return The same library but overwritten with the new values. A dataframe containing all details on the removed and added ids.
#'
#' @export

blib_seamlessIDswap <- function(db_main, db_redundant, rt, mz, tol= 1e-12, file= NA_character_){

  #check required packages
    if (!requireNamespace("DBI", quietly = TRUE)) {
      cli::cli_abort("Package 'DBI' is required but not installed. Please install it using BiocManager::install('DBI').")
    }
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    cli::cli_abort("Package 'RSQLite' is required but not installed. Please install it using BiocManager::install('RSQLite').")
  }


  # Check if the database file exists
  if (!file.exists(db_main)) { cli::cli_abort(c("x" = "The spectral library {.arg {db_main}} does not exist",
                                                      "i" = "Check the spelling or set the right working directory.")
                                                    )}

  if (!file.exists(db_redundant)) {cli::cli_abort(c("x" = "The spectral library {.arg{db_redundant}} does not exist",
                                                      "i" = "Check the spelling or set the right working directory.")
                                                      )}

  #databases are super-assigned to be used by all underlying functions.
  conn_redun   <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_redundant)
  conn_main    <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_main)

  #remove databases from global environment
  on.exit({
    DBI::dbDisconnect(conn_redun)
    DBI::dbDisconnect(conn_main)
  }, add = TRUE)



  DBI::dbBegin(conn_main)
  DBI::dbBegin(conn_redun)


  tryCatch({


    df_list <-
      purrr::pmap(
        .l = list(
          rt = rt,
          mz = mz,
          tol = tol,
          file = file
        ),
        .f = function(rt, mz, tol, file) {
          .fetch_ID(
            rt = rt,
            mz = mz,
            tol = tol,
            file = file,
            conn_redun = conn_redun,
            conn_main = conn_main
          )
        },
        .progress = TRUE
      ) |>
      dplyr::bind_rows()

  DBI::dbCommit(conn_main)
  DBI::dbCommit(conn_redun)

  return(df_list)




  }, error = function(e) {
    # Rollback transaction if an error occurs
    DBI::dbRollback(conn_main)
    DBI::dbRollback(conn_redun)
    cli::cli_alert_danger("Error: {e$message}")
    cli::cli_alert_info("ID swapping failed.")
    stop(e)
  })

}


#' @noRd
#top function
.fetch_ID <- function(rt, mz, tol, file,
                      conn_redun = conn_redun,
                      conn_main = conn_main){


  #get id for specific rt and mz from reduandant library
  id_table <- .blib_getRefID(conn= conn_redun, rt, mz, tol= 1e-12, file=file)


  id_redun  <- id_table$id
  id_mz     <- id_table$precursorMZ
  id_z      <- id_table$precursorCharge
  id_rt     <- trunc(100 * id_table$retentionTime) / 100
  id_modseq <- id_table$peptideModSeq
  id_file   <- id_table$fileName


  #transfer the id from the redundant library to the main library
    # Transfer ID data like rt, mz, ...
  .blib_transferID(conn_from = conn_redun, conn_to = conn_main, table = "RefSpectra", ids = id_redun)

    # Transfer intensity values of fragment ions
  .blib_transferID(conn_from = conn_redun, conn_to = conn_main, table = "RefSpectraPeaks", ids = id_redun )

  #make the id from the reduandant library as main id. i.e. best spectrum and get id_ref
  .blib_swapID(conn = conn_main, ids = id_redun) -> id_ref

  #remove initial id (that was representing the best sepctrum) from the main library
  .blib_removeID(conn = conn_main, ids = id_ref)


  df_summary <- dplyr::tbl(conn_redun, "RefSpectra") |>
    dplyr::collect() |>
    dplyr::filter( id %in% c(id_redun, id_ref)) |>
    dplyr::mutate(status = ifelse(id == id_redun, "added", "removed"), .before = 2)


  return(df_summary)



}


#' @noRd
#get file names
.blib_getFileID <- function(conn){

  src_file <- dplyr::tbl(conn, "SpectrumSourceFiles") |>
    dplyr::collect() |>
    dplyr::mutate(fileName = sapply(fileName, function(x) tools::file_path_sans_ext(basename(x)))) |>
    dplyr::select(-cutoffScore)

  return(src_file)

}


#' @noRd
#get id of an ID from library by providing mz and RT
.blib_getRefID <- function(conn, rt, mz, tol = 1e-12, file = NA_character_) {

    ID <- dplyr::tbl(conn, "RefSpectra") |>
    dplyr::collect() |>
    dplyr::filter(
      round(retentionTime, 1L) == round(rt, 1L) &
        dplyr::near(trunc(100 * precursorMZ) / 100, trunc(100 * mz) / 100, tol = tol)
    ) |>
    dplyr::select(
      id,
      precursorMZ,
      precursorCharge,
      peptideModSeq,
      retentionTime,
      fileID
    ) |>
    dplyr::left_join(.blib_getFileID(conn = conn), by = dplyr::join_by("fileID" == "id"))


if (!is.na(file) & nrow(ID) > 1L) {
    # Filter further if a specific file is provided
  ID <- ID |>
      dplyr::filter(fileName == file)

}else if (nrow(ID) == 0L) {
    cli::cli_abort("No matching ID was found for rt = {rt} and mz = {mz}.")

  }else if (is.na(file) & nrow(ID) > 1L) {
    cli::cli_abort(c(
      "The {.arg rt} and {.arg mz} must correspond to a unique ID.",
      "i" = "Specify {.arg file} name to disambiguate.",
      "x" = "Multiple IDs found having rt = {rt} and mz = {mz} in the redundant spectral library."
    ))
  }

    return(ID)

}





#' @noRd
#transfer ID from redundant library to main library
.blib_transferID <- function(conn_from, conn_to, table = c("RefSpectra", "RefSpectraPeaks"), ids) {


  id_col <- ifelse(table == "RefSpectra", "id", "RefSpectraID")

  if(id_col == "id"){
    # extract row
    rows_to_copy <- dplyr::tbl(conn_from, table) |>
      dplyr::filter(!!rlang::sym(id_col) %in% ids) |>
      dplyr::collect()

    #extract (modified)sequence & charge to group_by
    modseq <- rows_to_copy |>
      dplyr::pull(peptideModSeq)

    charge <- rows_to_copy |>
      dplyr::pull(precursorCharge)

    #modify the copy number in 'copies' column. For seamless swap of bestSpectrum ID.
    copy_update <- dplyr::tbl(conn_from, table) |>
      dplyr::filter(peptideModSeq == modseq, precursorCharge == charge) |>
      dplyr::summarise(sum_copies = sum(copies), .by = c('peptideModSeq', 'precursorCharge')) |>
      dplyr::pull(sum_copies)

    rows_to_copy <- rows_to_copy |>
      dplyr::mutate(copies = copy_update)

  }else{
    #no need to do anything when transfer is from 'RefSpectraPeaks' table
    rows_to_copy <- dplyr::tbl(conn_from, table) |>
      dplyr::filter(!!rlang::sym(id_col) %in% ids) |>
      dplyr::collect()
  }



  #Check if any rows exist
  if (nrow(rows_to_copy) == 0) {

   cli::cli_abort("ID was not found in the {.arg conn_from}.")

  }else{

    DBI::dbWriteTable(conn = conn_to, name = table, value = rows_to_copy, overwrite = FALSE, append = TRUE)

  }

}



#' @noRd
#This function operates on the retention time table of the main library to promote redun_id and silence ref_id
.blib_swapID <- function(conn, ids){

  rt_table <-  dplyr::tbl(conn, 'RetentionTimes') |>
    dplyr::collect()

  ref_id <- rt_table |>
    dplyr::filter(RedundantRefSpectraID == ids) |>
    dplyr::pull(RefSpectraID)

  rt_table_updated <-  rt_table |>
    dplyr::mutate(RefSpectraID = ifelse(
      RefSpectraID == ref_id,
      as.integer(ids),
      RefSpectraID
    ),
    bestSpectrum = dplyr::case_when(
      RedundantRefSpectraID == ids ~ 1L,
      RedundantRefSpectraID == ref_id ~ 0L,
      .default = as.integer(bestSpectrum)
    ))




  DBI::dbWriteTable(conn = conn, name= 'RetentionTimes', value = rt_table_updated, overwrite = TRUE)

  return(as.integer(ref_id))

}


#' @noRd
#remove id from main library RefSpectra and RefSpectraPeaks tables.
.blib_removeID <- function(conn, ids) {

  DBI::dbExecute(conn = conn, paste0("DELETE FROM RefSpectra WHERE id = ", ids))
  DBI::dbExecute(conn = conn, paste0("DELETE FROM RefSpectraPeaks WHERE RefSpectraID = ", ids))

}




