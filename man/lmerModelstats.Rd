% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmerModelstats.R
\name{lmerModelstats}
\alias{lmerModelstats}
\title{lmerModelstats}
\usage{
lmerModelstats(
  fit,
  altfit = NULL,
  df = NULL,
  confidenceinterval = 0.95,
  studywiseAlpha = 0.05
)
}
\arguments{
\item{fit}{lmer object.}

\item{altfit}{lmer object to use as comparison in hierarchical regression. If NULL then standard contrast against an empty model with just random effects is performed.}

\item{df}{text parameter to indicate what degrees of freedom to use. Options are Kenward-Roger (default) or Shattertwaite}

\item{confidenceinterval}{Decimal representation of confidence interval. Default 0.95.}

\item{studywiseAlpha}{Decimal representation of alpha level. Default 0.05.}
}
\value{
A list of standard reporting statistics
}
\description{
Obtain standard linear model reporting statistics for a lmer test
}
\examples{

    fit <- lmerTest::lmer(Alertness ~ Group + (1 | PartID), data = Rmimic::alertness)
    results <- Rmimic::lmerModelstats(fit)
    
    altfit <- lmerTest::lmer(Alertness ~ Time + (1 | PartID), data = Rmimic::alertness)
    fit <- lmerTest::lmer(Alertness ~ Group + Time + (1 | PartID), data = Rmimic::alertness)
    results <- Rmimic::lmerModelstats(fit, altfit)

}
\author{
Matthew B. Pontifex, \email{pontifex@msu.edu}, June 3, 2025
}
