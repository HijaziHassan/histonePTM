% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_scheduledPRMranges.R
\name{plot_schdueldPRMranges}
\alias{plot_schdueldPRMranges}
\title{Plot retention time window ranges for a PRM scheduled experiment.}
\usage{
plot_schdueldPRMranges(
  df,
  tr_start,
  tr_end,
  y_axis,
  label,
  save_plot = FALSE,
  ...
)
}
\arguments{
\item{df}{A dataframe containing at least three columns: m/z and start and end retention times.}

\item{tr_start}{starting retention time column name.}

\item{tr_end}{ending retention time column name.}

\item{y_axis}{column name of the variable to appear on x-axis. Default is "m/z" but the axis texts are balnked out.}

\item{label}{Default is m/z. number will appear in the middle of the horizontal bar.}

\item{save_plot}{\code{logical} in case needed to save the plot}

\item{...}{For any additions}
}
\value{
plot with ranges representing intervals when each precursor m/z is monitored.
}
\description{
Plot retention time window ranges for a PRM scheduled experiment.
}
