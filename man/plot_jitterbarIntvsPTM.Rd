% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_jitterbarIntvsPTM.R
\name{plot_jitterbarIntvsPTM}
\alias{plot_jitterbarIntvsPTM}
\title{Intensity vs PTM grouped bar plot with jittered points}
\usage{
plot_jitterbarIntvsPTM(
  dataset,
  x_axis,
  y_axis,
  condition,
  id_col,
  plot_title = NULL,
  max_cutoff = NULL,
  fun = c("mean", "median"),
  error_type = c("none", "CI", "SE", "SD"),
  conf_level = 0.95,
  scale = 100,
  cond_order = NULL,
  save_plot = FALSE,
  output_dir = NULL
)
}
\arguments{
\item{dataset}{A dataframe in long format with at least 4 columns: Intensity, PTM, sequence (or sequence label) and Condition. check the arguments below.
A 5th column can be provided to be the plot title.}

\item{x_axis}{x variable (PTM column)}

\item{y_axis}{y variable (intensity column). If already values are percentage, set \code{scale} to 1.}

\item{condition}{The condition column (WT vs disease, concentration, ...)}

\item{id_col}{unique ID column such as sequence or sequence label}

\item{plot_title}{(\code{optional})A column with values to be the plot_title (optional).}

\item{max_cutoff}{(\code{optional; numeric} A the maximum value above which values will be filtered out. It is added to draw separately low abundant PTMs.}

\item{fun}{mean (\code{default}) or median. This will be the height of the bar.}

\item{error_type}{one of "CI" (confidence interval), "SE" (standard error), "SD" (standard deviation).}

\item{conf_level}{Confidence level (e.g. 0.95 (default), 0.99, etc...)}

\item{scale}{100 (\code{default}). If you want to keep values as they are use 1. No other values are allowed.}

\item{cond_order}{(optional). A character vector containing conditions according to which bars will be ordered in the plot.}

\item{save_plot}{(\code{logical}; FALSE \code{default})}

\item{output_dir}{(character; Optional) The output directory for the plots to be saved in. Default to working directory if unassigned.}
}
\value{
A dataframe grouped by `id_col` (and `plot_title` if passed) with a nested list-column harboring the generated bar plots
(representing either `mean` or `median`) with jitter points (corresponding to individual measurements per `condition`).
}
\description{
Individual measurements will be plotted as data points. Their median (or mean) will be represented as a bar. Each bar will be colored according to the condition.
}
