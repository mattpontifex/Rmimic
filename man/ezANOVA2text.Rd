% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ezANOVA2text.R
\name{ezANOVA2text}
\alias{ezANOVA2text}
\title{ezANOVA2text}
\usage{
ezANOVA2text(
  result,
  numparticipants = NULL,
  feffect = NULL,
  sphericity = NULL,
  confidenceinterval = 0.95,
  studywiseAlpha = 0.05
)
}
\arguments{
\item{result}{ezANOVA output}

\item{numparticipants}{Number of unique participants in ANOVA}

\item{feffect}{Parameter to select effect size estimate. Default is Partial Eta Squared. Other option is Generalized Eta Squared.}

\item{sphericity}{Parameter to select Sphericity correction method. Default is Greenhouse-Geisser. Other option is Huynh-Feldt.}

\item{confidenceinterval}{Decimal representation of confidence interval. Default 0.95.}

\item{studywiseAlpha}{Decimal representation of alpha level. Default 0.05.}
}
\value{
A list with the same elements as the ezANOVA output:
\item{aov}{ANOVA object}
\item{ANOVA}{ANOVA table}
}
\description{
Output ezANOVA results in APA style format with effect sizes and confidence intervals.
}
\author{
Matthew B. Pontifex, \email{pontifex@msu.edu}, October 29, 2019
}
