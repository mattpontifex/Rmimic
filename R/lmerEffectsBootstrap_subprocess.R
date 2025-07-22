#' lmerEffectsBootstrap_subprocess
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerEffects function. 
#'
#' @param results List containing output from lmerEffects
#' @param repetitions Numeric entry indicating the number of repetitions to perform
#' @param method String entry indicating the approach. Default simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. This keeps the data around the original mean and standard deviation, thus as sample size increases the effect size will also increase. Resample performs data resampling with replacement from the existing dataset. Parametric simulates data from a multivariate normal distribution using the MASS::mvrnorm function. This keeps the effect size the same as it allows the group variation to increase with larger samples. Nonparametric simulates data from a multivariate nonnormal distribution using the mnonr::unonr function. This keeps the effect size the same as it allows the group variation to increase with larger samples.
#' @param subsample Numeric entry 0 to 1 indicating what percentage of the data to use for informing the parametric simulation. Ignored for method resample.
#' @param inflation Numeric entry indicating how many times the original sample to simulate too. Ignored for method resample.
#' @param resample_min Numeric entry indicating the minimum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param resample_max Numeric entry indicating the maximum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' 
#' @return
#' \item{results}{List of outputs.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom dplyr slice_sample
#' @importFrom progressor with_progress
#' @importFrom stats runif model.frame
#'
#' @export

lmerEffectsBootstrap_subprocess <- function(results, repetitions, resample_min, resample_max, subsample, inflation, method, boolposthoc, p) {
  
  # Define function for a single repetition
  run_one <- function(results, resample_min, resample_max, subsample, inflation, method, boolposthoc) {
    # populates dataset subsampled from the original sample with replacement
    if (method == "resample") {
      # Determine how large a random sample
      numberofsamples <- floor(stats::runif(1, min=resample_min, max=resample_max))
      smp <- dplyr::slice_sample(stats::model.frame(results$fit), n=numberofsamples, replace=TRUE)
    }
    if (method == "parametric") {
      # Useful for growing the sample 
      # keeps the exact same effect size by allowing the group variation to increase with larger samples
      smp <- Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=TRUE, method = "covariance")
    }
    if (method == "nonparametric") {
      # Useful for growing the sample 
      # keeps the exact same effect size by allowing the group variation to increase with larger samples
      smp <- Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=FALSE, method = "covariance")
    }
    if (method == "default") {
      # works the best as it is a wrapper around simulate - not ideal for growing the sample as the effect size will grow with it but will keep the data around the original mean and standard deviation
      smp <- Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, method = "conditionaldistribution")
    }
    
    # rerun model on new data
    newfit <- tryCatch({
      newfit <- invisible(suppressWarnings(suppressMessages(update(results$fit, data=smp, evaluate = TRUE))))
    }, error = function(e) {
      cat(sprintf('lmerEffectsBootstrap - model failure\n'))
      newfit <- NULL
    })
    
    
    if (!is.null(newfit)) {
      # extract information
      newresults <- invisible(suppressWarnings(suppressMessages(lmerEffects(newfit, dependentvariable=results$dependentvariable, subjectid=results$subjectid, within=results$within, df = results$df, confidenceinterval=results$confidenceinterval, studywiseAlpha=results$studywiseAlpha, suppresstext=TRUE, smp=smp))))
      # compute posthoc if previously run - function should have stored the necessary information
      if (boolposthoc) {
        newresults <- invisible(suppressWarnings(suppressMessages(lmerPosthoc(newresults, between=results$between, within=results$within, covariates=results$covariates, planned=results$stats$Effect, posthoclimit=results$posthoclimit, calltype='subprocess', posthoccorrection='none'))))
      }
      newresults$descriptives <- lmerEffects_simpledesc(newresults) # store these
      newresults$fit <- NULL
      return(newresults)
    } else {
      return(NULL)
    }
  }
  
  if (!is.null(method)) {
    if (toupper(method) == toupper("resample")) {
      method = "resample"
    } else if (toupper(method) == toupper("parametric")) {
      method = "parametric"
    } else if (toupper(method) == toupper("nonparametric")) {
      method = "nonparametric"
    } else if (toupper(method) == toupper("default")) {
      method = "default"
    } else {
      method = "default"
    }
  } else {
    method = "default"
  }
  
  library(future)
  library(future.apply)
  handlers(global = TRUE)
  
  # Set up parallel plan
  plan(multisession, workers = availableCores() - 1)
  
  resstore <- list()
  
  # Wrap the call in `with_progress()`
  resstore <- with_progress({
    future_lapply(1:repetitions, function(i) {
      result <- run_one(results, resample_min, resample_max, subsample, inflation, method, boolposthoc)
      p()
      result
    }, future.seed = TRUE)
  })
  
  return(resstore)
}

