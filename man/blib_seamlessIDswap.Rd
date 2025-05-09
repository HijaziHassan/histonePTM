% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blib_seamlessIDswap.R
\name{blib_seamlessIDswap}
\alias{blib_seamlessIDswap}
\title{Change best spectrum ID in \emph{.blib} spectral library}
\usage{
blib_seamlessIDswap(
  db_main,
  db_redundant,
  rt,
  mz,
  tol = 1e-12,
  file = NA_character_
)
}
\arguments{
\item{db_main}{The name of the main Skyline \code{.blib} library file (\code{string}).}

\item{db_redundant}{The name of the redundant Skyline \code{.blib} library file (\code{string}).}

\item{rt}{Retention time (\code{numeric}) of the ID to be set as the best spectrum (rounded to first decimal).
Also accepts numeric vector.}

\item{mz}{\emph{m/z} (\code{numeric}) value of the target ID (it will be truncated to two decimal places).
Also accepts numeric vector.}

\item{tol}{The tolerance (\code{optional}) to be used to match the provided \emph{m/z} (\code{numeric}).
The default value is used to have a perfect match.}

\item{file}{MS raw file name without extension (\code{optional}). In case of multiple matches to \code{rt},
 it is necessary to provide the \code{file} name to force specific ID selection (the file name is found at the top of
chromatogram window of each sample in Skyline interface).
Also accepts character vector.}
}
\value{
The same library but overwritten with the new values. A dataframe containing all details on the removed and added ids.
}
\description{
To swap IDs between redundant library and main library in case of wrong peak picking or (near)-coeluting isobaric peptides
in Skyline software.
}
