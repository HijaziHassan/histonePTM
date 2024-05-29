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
