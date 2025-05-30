% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_intDensity.R
\name{plot_intDensity}
\alias{plot_intDensity}
\title{A ridgeline visualization of intensities per sample}
\usage{
plot_intDensity(
  df,
  seq_col,
  ptm_col,
  int_col,
  format = c("wide", "long"),
  save_plot = FALSE
)
}
\arguments{
\item{df}{Dataset in wide format.}

\item{seq_col}{sequence column.}

\item{ptm_col}{modification column.}

\item{int_col}{The intensity columns. Either a vector `c('column1', 'column2', ...)` or using tidyselect functions (e.g. starts_With('abundance'))}

\item{format}{'wide' or 'long'.}

\item{save_plot}{\code{bool}, FALSE (default).}
}
\value{
ridgeline plot
}
\description{
A plot the gives an eye-bird view of the intensity distribution in each sample.
}
