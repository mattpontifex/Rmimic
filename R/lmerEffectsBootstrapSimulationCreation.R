#' lmerEffectsBootstrapSimulationCreation
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerEffects function. This portion handles model creation as large repetitions may take substantial time. This function is additive so it can allow the process to be performed in chunks if needed. 
#'
#' @param results List containing output from lmerEffects
#' @param repetitions Numeric entry indicating the number of repetitions to perform
#' @param method String entry indicating the approach. Default simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. This keeps the data around the original mean and standard deviation, thus as sample size increases the effect size will also increase. Resample performs data resampling with replacement from the existing dataset. Parametric simulates data from a multivariate normal distribution using the MASS::mvrnorm function. This keeps the effect size the same as it allows the group variation to increase with larger samples. Nonparametric simulates data from a multivariate nonnormal distribution using the mnonr::unonr function. This keeps the effect size the same as it allows the group variation to increase with larger samples.
#' @param subsample Numeric entry 0 to 1 indicating what percentage of the data to use for informing the parametric simulation. Ignored for method resample.
#' @param inflation Numeric entry indicating how many times the original sample to simulate too. Ignored for method resample.
#' @param resample_min Numeric entry indicating the minimum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param resample_max Numeric entry indicating the maximum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param tmpdir Path to write temporary model files to
#'
#' @return
#' \item{results}{List containing the original output from lmerEffects with some additional tag information.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 23, 2025
#'
#' @importFrom stats model.frame runif
#' @importFrom progressor progressor with_progress
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom dplyr slice_sample
#' 
#' @examples
#' \dontrun{
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
#'                covariates=NULL, planned=c('Group'), posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
#'     results <- Rmimic::lmerEffectsBootstrapSimulationCreation(results, repetitions=999, tmpdir='/tempdirectory/')
#'     }
#'
#' @export

lmerEffectsBootstrapSimulationCreation <- function(results, repetitions, resample_min=NULL, resample_max=NULL, subsample=0.96, inflation=1.0, method='default', tmpdir=NULL) {
  
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
      newresults$fit <- NULL # clear model fit to save space
      return(newresults)
    } else {
      return(NULL)
    }
  }
  
  
  
  
  
  library(future)
  library(future.apply)
  library(progressr)
  
  
  
  totalsample <- nrow(stats::model.frame(results$fit))
  if (is.null(resample_min)) {
    resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (is.null(resample_max)) {
    resample_max <- nrow(stats::model.frame(results$fit))
  }  
  methodofposthoccorrection <- 'false'
  if ('posthoccorrection' %in% names(results)) {
    methodofposthoccorrection <- results$posthoccorrection
  }
  boolposthoc <- FALSE
  if ('posthoc' %in% names(results)) {
    boolposthoc <- TRUE
    # rerun but force all posthocs
    results <- invisible(suppressWarnings(suppressMessages(lmerPosthoc(results, between=results$between, within=results$within, covariates=results$covariates, planned=results$stats$Effect, posthoclimit=results$posthoclimit, posthoccorrection='none', progressbar=FALSE))))
  }
  results$descriptives <- lmerEffects_simpledesc(results) # store these
  
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
  
  
  # Set up parallel plan
  plan(multisession, workers = availableCores() - 1)
  
  # Parameters
  max_files <- repetitions
  repetitions <- max_files*1.25  # ask for more
  dir.create(tmpdir, showWarnings = FALSE)
  
  # Check how many result files exist
  file_list <- list.files(tmpdir, pattern = "^result_.*\\.RData$")
  # Update to reflect if any files have already been created
  repetitions <- repetitions - length(file_list) + 1
  
  # Enable progress handler
  #handlers("txtprogressbar") 
  handlers(list(
    handler_progress(
      format   = "lmerEffectsBootstrapSimulationCreation() [:bar] :percent :eta",
      width    = 120,
      complete = "="
    )
  ))
  
  # Loop and track
  captureout <- with_progress({
    p <- progressor(along = 1:repetitions)
    
    future_lapply(1:repetitions, function(i) {
      # Check how many result files exist
      file_list <- list.files(tmpdir, pattern = "^result_.*\\.RData$")
      
      if (length(file_list) < max_files) {
        result <- run_one(results, resample_min, resample_max, subsample, inflation, method, boolposthoc)
        if (!is.null(result)) {
          save(result, file = tempfile(pattern = "result_", tmpdir = tmpdir, fileext = ".RData"))
        }
      }
      p()
    }, future.seed = TRUE)
  })
  
  # check final result
  Sys.sleep(1) # to make sure files are written
  
  file_list <- list.files(tmpdir, pattern = "^result_.*\\.RData$")
  if (length(file_list) > max_files) {
    # created too many files
    numtoremove <- length(file_list) - max_files
    files_to_delete <- file.path(tmpdir, sample(file_list, numtoremove, replace=FALSE))
    unlink(files_to_delete)
  }
  
  # include some tags to pass along information
  results$futuretag <- list()
  results$futuretag$repetitions <- max_files
  results$futuretag$resample_min <- resample_min
  results$futuretag$resample_max <- resample_max
  results$futuretag$subsample <- subsample
  results$futuretag$inflation <- inflation
  results$futuretag$method <- method
  results$futuretag$tmpdir <- tmpdir
  results$futuretag$boolposthoc <- boolposthoc
  results$futuretag$methodofposthoccorrection <- methodofposthoccorrection
  results$futuretag$totalsample <- totalsample
  
  Sys.sleep(1) # to make sure files are written
  
  return(results)
}
