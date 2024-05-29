



#' Annotate MS/MS Spectra
#' @description
#' wrapper function around functions from \code{spectrum_utils} library in Python.
#'
#' @param Profrma_peptide Peptide in proforma representation
#' @param prec_mz m/z of the precursor ion
#' @param prec_z charge state of the precursor ion
#' @param mz \code{numeric} vector containing mz values
#' @param intensity \code{numeric} vector containing intensity values
#' @param title title of the annotate MS/MS plot
#' @param output_plot_name Name of the plot in case you want to save it
#'
#'
#' @return Annotate MS/MS plot
#' @export
plot_annotateSpectrum <- function(Profrma_peptide,
                                  prec_mz, prec_z,
                                  mz, intensity,
                                  title,
                                  output_plot_name){


reticulate::source_python("./inst/python/plot_annotateSpectra.py")


plot_annotateSpectra(identifier= title,
                       peptide= Profrma_peptide,
                       precursor_mz= prec_mz,
                       precursor_charge = prec_z,
                       mz = mz,
                       intensity = intensity,
                       output_plot_name = output_plot_name)
}


