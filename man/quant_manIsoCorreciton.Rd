% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant_manIsoCorreciton.R
\name{quant_manIsoCorreciton}
\alias{quant_manIsoCorreciton}
\title{Correct of Naturally Occuring Abundances}
\usage{
quant_manIsoCorreciton(df, time_col, t0, cond_col, avg_by)
}
\arguments{
\item{df}{A dataframe as Picor output but with raw values.}

\item{time_col}{Name of thetime column.}

\item{t0}{Value of the t0 ( 0 , t0, T0 ...)}

\item{cond_col}{Condition column (Concentration, disease vs non-disease ...)}
}
\value{
A dataframe with corrected and normalized values.
}
\description{
Subtract naturally occuring abundances of MS XICs in stable isotope tracing experiments
}
