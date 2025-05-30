% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmerPosthoc.R
\name{lmerPosthoc}
\alias{lmerPosthoc}
\title{lmerPosthoc}
\usage{
lmerPosthoc(
  results,
  between = NULL,
  within = NULL,
  covariates = NULL,
  dependentvariable = NULL,
  subjectid = NULL,
  df = NULL,
  planned = NULL,
  suppress = NULL,
  posthoccorrection = NULL,
  confidenceinterval = NULL,
  studywiseAlpha = NULL,
  posthoclimit = 6,
  calltype = NULL,
  progressbar = TRUE
)
}
\arguments{
\item{results}{list output from the lmerEffects function}

\item{between}{text list specifying the between subjects variables}

\item{within}{text list specifying any variables that should be considered as within subjects. If null assumes all variables are between subjects.}

\item{covariates}{text list specifying the variables that should be considered as non interactive covariates.}

\item{planned}{text list specifying any effect to show the post-hoc comparisons even if they are not significant.}

\item{suppress}{text list specifying any effect to ignore even if it is significant.}

\item{posthoccorrection}{text parameter to indicate what post-hoc comparisons should be performed. Default is False Discovery Rate Control. Other options are Bonferroni, Holm-Bonferroni, Sidak. None will skip post hoc corrections.}

\item{confidenceinterval}{Decimal representation of confidence interval. Default 0.95.}

\item{studywiseAlpha}{Decimal representation of alpha level. Default 0.05.}

\item{posthoclimit}{integer specifying how many factors in an interaction should be broken down. Default is 6 indicating posthoc results will be provided for up to a 6 way interaction.}

\item{progressbar}{boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.}
}
\value{
results list output from the lmerEffects function
\item{posthoc}{list output containing the posthoc test results}
}
\description{
Posthoc decomposition of univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. Interactions are decomposed multiple ways (A holding B, B holding A) and superseeding interactions suppress lower level effects tests (no posthoc test of A:B if A:B:C is significant). Tests of A:B can still be obtained using the planned parameter if desired.
}
\examples{

    fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
    results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
    results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
               covariates=NULL, planned=c('Group'), posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
               
}
\author{
Matthew B. Pontifex, \email{pontifex@msu.edu}, May 1, 2025
}
