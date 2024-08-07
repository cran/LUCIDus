% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.lucid_serial}
\alias{summary.lucid_serial}
\title{Summarize results of the serial LUCID model}
\usage{
\method{summary}{lucid_serial}(object, ...)
}
\arguments{
\item{object}{A LUCID model fitted by \code{\link{estimate_lucid}}}

\item{...}{Additional argument of \code{boot.se}, which is an object returned by \code{\link{boot_lucid}}}
}
\value{
A list, containing the extracted key parameters from the LUCID model that can be used to print the summary table
}
\description{
Summarize results of the serial LUCID model
}
