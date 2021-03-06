% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table2console.R
\name{table2console}
\alias{table2console}
\title{table2console}
\usage{
table2console(
  dataframein,
  sepgap = NULL,
  spansize = 95,
  headers = TRUE,
  alternate = FALSE,
  seperators = FALSE
)
}
\arguments{
\item{dataframein}{Data frame to be printed}

\item{sepgap}{Data frame of the same size and with the same headers as the dataframein specifying the number of characters for each column.}

\item{spansize}{Number of characters on a line}

\item{headers}{Parameter to control if the variable labels should be printed. Default is TRUE.}

\item{alternate}{Parameter to control if the first column should have its own row. Default is FALSE.}

\item{seperators}{Parameter to control if the seperator lines should be shown. Default is FALSE.}
}
\description{
Function that takes a data table and prints a nicely formatted table to the console
}
\author{
Matthew B. Pontifex, \email{pontifex@msu.edu}, October 29, 2019
}
