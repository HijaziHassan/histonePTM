
#' Annotate MS/MS Spectra
#' @description
#' wrapper function around functions from \code{spectrum_utils} library in Python.
#'
#' @param Profrma_peptide Peptide in proforma representation
#' @param prec_mz m/z of the precursor ion (double)
#' @param prec_z charge state of the precursor ion (integer)
#' @param mz \code{numeric} vector containing mz values
#' @param intensity \code{numeric} vector containing intensity values
#' @param annot_diag (bool, FALSE (default)) Annotate CyIm ions produced from modified lysines.
#' @param nloss (bool, TRUE (default)). Annotate neutral losses (-NH3 by "*" and -H2O by "o")
#' @param min_int (double, 0.01 (default)). To filter intensity which are `min_int` of the base peak.
#' @param ion_types Fragment ion types (i.e. "by' (Default))
#' @param tol_mz fragment mass tolerance (e.g, 20 (Default))
#' @param tol_mode 'ppm' (Default) or 'Da'
#' @param title title of the annotate MS/MS plot
#' @param output_plot_name Name of the plot in case you want to save it
#'
#' @return Annotated MS/MS plot
#' @export
plot_annotateSpectrum <- function(Profrma_peptide,
                                  prec_mz,
                                  prec_z,
                                  mz, intensity,
                                  title= "",
                                  annot_diag= FALSE,
                                  nloss= TRUE,
                                  min_int= 0.01,
                                  ion_types= "by",
                                  tol_mz = "20",
                                  tol_mode = "ppm",
                                  output_plot_name= ""){


  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  My_Package_Name_path <- system.file("python", package = "histonePTM")
  reticulate::py_run_string(paste0("import sys; sys.path.append('", My_Package_Name_path, "')"))

  reticulate::source_python(system.file("python",
                                        "plot_annotateSpectra.py",
                                        package = "histonePTM",
                                        mustWork = TRUE
  ))


  reticulate::py$plot_annotateSpectra(
    identifier = title,
    peptide = Profrma_peptide,
    precursor_mz = prec_mz,
    precursor_charge = prec_z,
    mz = mz,
    intensity = intensity,
    nloss = nloss,
    annot_diag = annot_diag,
    ion_types = ion_types,
    min_int = min_int,
    tol_mz = tol_mz,
    tol_mode = tol_mode,
    output_plot_name = output_plot_name
  )
}


