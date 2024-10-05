#' Create Folder
#'
#' @description
#' Create a folder carrying the name of the analyzed file (.csv or .xlsx). If a path is provided, only what after the last slash is considered.
#'
#' @param foldername name of the folder with or without forward-slash (/). Only the name after the last slash will be used. Any extension will be removed.
#' @importFrom cli cli_alert_info
#' @return folder
#' @export
misc_createFolder <- function(foldername){

  foldername= tools::file_path_sans_ext({{foldername}})
  if(!dir.exists(file.path(foldername))){
    dir.create(file.path(foldername))
    foldername = sub(".*/(.+)$", "\\1", foldername)
    cli::cli_inform("'{foldername}' folder is created.")

  }else{
    foldername = sub(".*/(.+)$", "\\1", foldername)
    cli::cli_alert_info("\nHey, '{foldername}' folder already exists!")}
}

#' @examples
#' misc_createFolder('my_folder')
