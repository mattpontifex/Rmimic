% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell2span.R
\name{cell2span}
\alias{cell2span}
\title{cell2span}
\usage{
cell2span(dataframein, sepgap = NULL, spansize = 95, trimwarn = TRUE)
}
\arguments{
\item{dataframein}{Data frame to be transformed}

\item{sepgap}{Data frame of the same size and with the same headers as the dataframein specifying the number of characters for each column.}

\item{spansize}{Number of characters on a line}

\item{trimwarn}{Parameter to control if user should be warned that text was trimmed to fit. Default is TRUE}
}
\description{
Function that takes a vector of data table cells and creates a character string of a particular length. The text empty will result in an open span in the output
}
\author{
Matthew B. Pontifex, \email{pontifex@msu.edu}, October 29, 2019
}
