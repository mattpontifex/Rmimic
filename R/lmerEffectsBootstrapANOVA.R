#' lmerEffectsBootstrapANOVA
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerPosthoc function. This portion handles the primary anova tables.
#'
#' @param results List containing output from lmerEffects
#' @param repetitions Numeric entry indicating the number of repetitions to perform
#' @param method String entry indicating the approach. Default simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. This keeps the data around the original mean and standard deviation, thus as sample size increases the effect size will also increase. Resample performs data resampling with replacement from the existing dataset. Parametric simulates data from a multivariate normal distribution using the MASS::mvrnorm function. This keeps the effect size the same as it allows the group variation to increase with larger samples. Nonparametric simulates data from a multivariate nonnormal distribution using the mnonr::unonr function. This keeps the effect size the same as it allows the group variation to increase with larger samples.
#' @param subsample Numeric entry 0 to 1 indicating what percentage of the data to use for informing the parametric simulation. Ignored for method resample.
#' @param inflation Numeric entry indicating how many times the original sample to simulate too. Ignored for method resample.
#' @param resample_min Numeric entry indicating the minimum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param resample_max Numeric entry indicating the maximum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param tmpdir Path to write temporary model files to
#' @param average text parameter to indicate how the results should be collapsed. Default is using the median. Other options is mean.
#' @param reporteddata text parameter to indicate if the posthoc text reports should use actual data (actual) or should report the simulated or resampled data means.
#'
#' @return
#' \item{results}{List containing bootstrapped output.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, August 18, 2025
#'
#' @importFrom stats model.frame runif
#' @importFrom progressr progressor with_progress
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom dplyr slice_sample
#' @importFrom data.table fwrite rbindlist fread setDT copy
#' @importFrom Rmisc CI
#' @importFrom miscTools colMedians
#' 
#' @export

lmerEffectsBootstrapANOVA <- function(results, ...) {
  options(warn = -1)
  
  bootstrap <- tryCatch({
    bootstrap <- list(...) 
  }, error = function(e) {
    bootstrap <- list() 
  })
  if (!('repetitions' %in% tolower(names(bootstrap)))) {
    bootstrap$repetitions <- 100
  }
  if (!('resample_min' %in% tolower(names(bootstrap)))) {
    bootstrap$resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (!('resample_max' %in% tolower(names(bootstrap)))) {
    bootstrap$resample_max <- nrow(stats::model.frame(results$fit))
  }
  if (!('subsample' %in% tolower(names(bootstrap)))) {
    bootstrap$subsample <- 1.0
  }
  if (!('inflation' %in% tolower(names(bootstrap)))) {
    bootstrap$inflation <- 1.0
  }
  if (!('method' %in% tolower(names(bootstrap)))) {
    bootstrap$method <- 'default'
  }
  if (!('tmpdir' %in% tolower(names(bootstrap)))) {
    bootstrap$tmpdir <- ''
  }
  if (!('average' %in% tolower(names(bootstrap)))) {
    bootstrap$average <- 'median'
  }
  if (!('reporteddata' %in% tolower(names(bootstrap)))) {
    bootstrap$reporteddata <- 'simulated'
  }
  
  # Define function for a single repetition
  run_one <- function(results, resample_min, resample_max, subsample, inflation, method) {
    # populates dataset subsampled from the original sample with replacement
    smp <- NULL
    if (method == "resample") {
      # Determine how large a random sample
      numberofsamples <- floor(stats::runif(1, min=resample_min, max=resample_max))
      smp <- invisible(suppressWarnings(suppressMessages(dplyr::slice_sample(stats::model.frame(results$fit), n=numberofsamples, replace=TRUE))))
    }
    if (method == "parametric") {
      # Useful for growing the sample 
      # keeps the exact same effect size by allowing the group variation to increase with larger samples
      smp <- invisible(suppressWarnings(suppressMessages(Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=TRUE, method = "covariance"))))
    }
    if (method == "nonparametric") {
      # Useful for growing the sample 
      # keeps the exact same effect size by allowing the group variation to increase with larger samples
      smp <- invisible(suppressWarnings(suppressMessages(Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=FALSE, method = "covariance"))))
    }
    if (method == "default") {
      # works the best as it is a wrapper around simulate - not ideal for growing the sample as the effect size will grow with it but will keep the data around the original mean and standard deviation
      smp <- invisible(suppressWarnings(suppressMessages(Rmimic::lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, method = "conditionaldistribution"))))
    }
    
    if (!is.null(smp)) {
      # rerun model on new data
      newfit <- tryCatch({
        newfit <- invisible(suppressWarnings(suppressMessages(update(results$fit, data=smp, evaluate = TRUE))))
      }, error = function(e) {
        cat(sprintf('lmerSimulateData - model failure\n'))
        newfit <- NULL
      })
    } else {
      newfit <- NULL
    }
    
    if (!is.null(newfit)) {
      # extract information
      newresults <- invisible(suppressWarnings(suppressMessages(lmerEffects(newfit, dependentvariable=results$dependentvariable, subjectid=results$subjectid, within=results$within, df = results$df, confidenceinterval=results$confidenceinterval, studywiseAlpha=results$studywiseAlpha, suppresstext=TRUE, smp=smp))))
      newresults$descriptives <- lmerEffects_simpledesc(newresults) # store these
      outresults <- list()
      outresults$descriptives <- newresults$descriptives
      outresults$stats <- newresults$stats
      outresults$randomstats <- newresults$randomstats
      outresults$rsquared <- newresults$rsquared
      return(outresults)
    } else {
      return(NULL)
    }
  }
  
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future)))))
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future.apply)))))
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(progressr)))))
  
  # tag model
  results$bootstrap <- list()
  results$bootstrap$repetitions <- bootstrap$repetitions
  results$bootstrap$resample_min <- bootstrap$resample_min
  results$bootstrap$resample_max <- bootstrap$resample_max
  results$bootstrap$subsample <- bootstrap$subsample
  results$bootstrap$inflation <- bootstrap$inflation
  results$bootstrap$method <- bootstrap$method
  results$bootstrap$tmpdir <- bootstrap$tmpdir
  results$bootstrap$average <- bootstrap$average
  results$bootstrap$reporteddata <- bootstrap$reporteddata
  
  bootstrap$totalsample <- nrow(stats::model.frame(results$fit))
  if (is.null(bootstrap$resample_min)) {
    bootstrap$resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (is.null(bootstrap$resample_max)) {
    bootstrap$resample_max <- nrow(stats::model.frame(results$fit))
  }  
  
  results$descriptives <- lmerEffects_simpledesc(results) # store these
  
  if (!is.null(bootstrap$method)) {
    if (toupper(bootstrap$method) == toupper("resample")) {
      bootstrap$method = "resample"
    } else if (toupper(bootstrap$method) == toupper("parametric")) {
      bootstrap$method = "parametric"
    } else if (toupper(bootstrap$method) == toupper("nonparametric")) {
      bootstrap$method = "nonparametric"
    } else if (toupper(bootstrap$method) == toupper("default")) {
      bootstrap$method = "default"
    } else {
      bootstrap$method = "default"
    }
  } else {
    bootstrap$method = "default"
  }
  
  # Set up parallel plan
  n_workers <- (availableCores() - 1)
  if (n_workers > 124) {
    n_workers <- 124
  }
  plan(multisession, workers = n_workers)
  
  # Parameters
  max_files <- bootstrap$repetitions
  repetitions <- floor(max_files+(n_workers*2))  # ask for more
  if (bootstrap$tmpdir == '') {
    #bootstrap$tmpdir <- file.path(getwd(), paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), collapse = ""))
    bootstrap$tmpdir <- tempfile(pattern = "anova", tmpdir = getwd(), fileext = "")
  }
  dir.create(bootstrap$tmpdir, showWarnings = FALSE)
  
  # Check how many result files exist
  file_list <- list.files(bootstrap$tmpdir, pattern = "^result_stats_.*\\.csv$")
  
  # Update to reflect if any files have already been created
  repetitions <- repetitions - length(file_list) + 1
  
  if (length(file_list) < max_files) {
    
    # Enable progress handler
    #handlers("txtprogressbar") 
    handlers(list(
      handler_progress(
        format   = "lmerEffectsBootstrapANOVA() simulating data [:bar] :percent :eta",
        width    = 120,
        complete = "="
      )
    ))
    
    # Loop and track
    captureout <- with_progress({
      step <- ceiling(repetitions / 10)
      p <- progressor(steps = 10)
      future_lapply(1:repetitions, function(i) {
        invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(mice)))))
        invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future)))))
        invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future.apply)))))
        invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(progressr)))))
        
        # Check how many result files exist
        file_list <- list.files(bootstrap$tmpdir, pattern = "^result_stats_.*\\.csv$")
        
        if (length(file_list) < max_files) {
          smp <- tryCatch({
            result <- run_one(results, bootstrap$resample_min, bootstrap$resample_max, bootstrap$subsample, bootstrap$inflation, bootstrap$method)
            if (!is.null(result)) {
              data.table::fwrite(result$descriptives, tempfile(pattern = "result_descriptives_", tmpdir = bootstrap$tmpdir, fileext = ".csv"), append = FALSE)
              data.table::fwrite(result$randomstats, tempfile(pattern = "result_randomstats_", tmpdir = bootstrap$tmpdir, fileext = ".csv"), append = FALSE)
              data.table::fwrite(result$rsquared, tempfile(pattern = "result_rsquared_", tmpdir = bootstrap$tmpdir, fileext = ".csv"), append = FALSE)
              data.table::fwrite(result$stats, tempfile(pattern = "result_stats_", tmpdir = bootstrap$tmpdir, fileext = ".csv"), append = FALSE)
            }
            smp <- NULL
          }, error = function(e) {
            smp <- NULL
          })
        }
        if ((i %% step == 0) || (i == repetitions)) {
          p()
        }
      }, future.seed = TRUE)
    })
    
  }
  
  # check final result
  Sys.sleep(1) # to make sure files are written
  
  # summarize descriptives
  file_list <- list.files(bootstrap$tmpdir, pattern = "^result_descriptives_.*\\.csv$", full.names=TRUE)
  if (length(file_list) > 100) {
    #datain <- data.table::rbindlist(future_lapply(file_list, data.table::fread), use.names=TRUE, fill=TRUE)
    # Enable progress handler
    #handlers("txtprogressbar") 
    handlers(list(
      handler_progress(
        format   = "lmerEffectsBootstrapANOVA() compiling data [:bar] :percent :eta",
        width    = 120,
        complete = "="
      )
    ))
    datain <- with_progress({
      step <- ceiling(length(file_list) / 10)        # 10% intervals
      p <- progressor(steps = 10)    # progress bar with 10 ticks
      data.table::rbindlist(
        future_lapply(seq_along(file_list), function(i) {
          out <- data.table::fread(file_list[i])
          if ((i %% step == 0) || (i == length(file_list))) {
            p()
          }
          out
        }),
        use.names = TRUE, fill = TRUE
      )
    })
  } else {
    datain <- data.table::rbindlist(lapply(file_list, data.table::fread), use.names=TRUE, fill=TRUE)
  }
  newstats <- results$descriptives
  for (cR in 1:nrow(newstats)) {
    tempdbs <- datain[Group == newstats$Group[cR]]
    colsofinterest <- c("N", "Missing", "Mean", "Median", "SD", "SE")
    for (cC in 1:length(colsofinterest)) {
      tempvect <- as.numeric(unlist(tempdbs[[colsofinterest[cC]]]))
      tempvect <- tempvect[which(!is.na(tempvect))]
      if (bootstrap$average == 'median') {
        newstats[cR,colsofinterest[cC]] <- median(tempvect, na.rm=TRUE)
      } else {
        newstats[cR,colsofinterest[cC]] <- mean(tempvect, na.rm=TRUE)
      }
    }
    
    distributiontable <- data.table::rbindlist(
      lapply(tempdbs$DistributionData, function(x) as.list(as.numeric(strsplit(x, ",")[[1]]))),
      fill = TRUE, use.names=FALSE
    )
    data.table::setDT(distributiontable)
    if (nrow(distributiontable) > 0) {
      if (bootstrap$average == "median") {
        newstats$DistributionData[cR] <- paste(miscTools::colMedians(as.matrix(distributiontable), na.rm=TRUE), collapse=",")
      } else {
        newstats$DistributionData[cR] <- paste(colMeans(distributiontable, na.rm=TRUE), collapse=",")
      }
    }
    
    # decision rule
    decisionindx <- tempdbs[DistributionDecision == "Normal", .N]
    if (decisionindx > 0) {
      if ((decisionindx / nrow(tempdbs)) >= 0.6) {
        newstats$DistributionDecision[cR] <- "Normal"
      } else {
        newstats$DistributionDecision[cR] <- "Not Normal"
      }
    }
  }
  
  # dirty but effective
  newstats$Mean <- round(round(round(round(as.numeric(newstats$Mean), digits=4), digits=3), digits=2), digits=1)
  newstats$Median <- round(round(round(round(as.numeric(newstats$Median), digits=4), digits=3), digits=2), digits=1)
  newstats$SD <- round(round(round(round(as.numeric(newstats$SD), digits=4), digits=3), digits=2), digits=1)
  newstats$SE <- round(round(round(round(as.numeric(newstats$SE), digits=4), digits=3), digits=2), digits=1)
  
  newstats$Mean <- sprintf('%.1f', newstats$Mean)
  newstats$Median <- sprintf('%.1f', newstats$Median)
  newstats$SD <- sprintf('%.1f', newstats$SD)
  newstats$SE <- sprintf('%.1f', newstats$SE)
  
  results$descriptives <- newstats
  rm(newstats)
  
  # summarize stats
  statsofinterst <- c('stats', 'randomstats', 'rsquared')
  for (cSi in 1:length(statsofinterst)) {
    
    file_list <- list.files(bootstrap$tmpdir, pattern = sprintf('^result_%s_.*\\.csv$', statsofinterst[cSi]), full.names=TRUE)
    if (length(file_list) > 100) {
      #datain <- data.table::rbindlist(future_lapply(file_list, data.table::fread), use.names=TRUE, fill=TRUE)
      handlers(list(
        handler_progress(
          format   = "lmerEffectsBootstrapANOVA() compiling data [:bar] :percent :eta",
          width    = 120,
          complete = "="
        )
      ))
      datain <- with_progress({
        step <- ceiling(length(file_list) / 10)        # 10% intervals
        p <- progressor(steps = 10)    # progress bar with 10 ticks
        data.table::rbindlist(
          future_lapply(seq_along(file_list), function(i) {
            out <- data.table::fread(file_list[i])
            if ((i %% step == 0) || (i == length(file_list))) {
              p()
            }
            out
          }),
          use.names = TRUE, fill = TRUE
        )
      })
    } else {
      datain <- data.table::rbindlist(lapply(file_list, data.table::fread), use.names=TRUE, fill=TRUE)
    }
    
    # copy to data.table
    textcall <- sprintf('newstats <- data.table::as.data.table(data.table::copy(results$%s))', statsofinterst[cSi])
    eval(parse(text=textcall))
    
    # clear to avoid confusion
    if ("textoutput" %in% names(newstats)) newstats[, textoutput := NA]
    if ("significance" %in% names(newstats)) newstats[, significance := NA]
    if ("F.value" %in% names(newstats)) newstats[, `:=`(F.value.ci.lower = NA, F.value.ci.upper = NA)]
    if ("p.value" %in% names(newstats)) newstats[, `:=`(p.value.ci.lower = NA, p.value.ci.upper = NA)]
    
    colsofinterest <- c("DFn", "DFd", "SSn", "SSd", "SSe", "F.value", "DF",
                        "LogLikelihood", "LRT", "p.value", "partialetasquared",
                        "fsquared","fsquared.ci.lower", "fsquared.ci.upper", "effects")
    
    for (cR in seq_len(nrow(newstats))) {
      
      if (all(c("SSn","SSd","F.value") %in% names(datain)) | all(c("LogLikelihood","LRT") %in% names(datain))) {
        tempdbs <- datain[Effect == newstats$Effect[cR]]
      } else if (all(c("portion","ci.lower","ci.upper") %in% names(datain))) {
        tempdbs <- datain[portion == newstats$portion[cR]]
      }
      
      # summarize numeric columns
      for (cC in colsofinterest) {
        if (cC %in% names(tempdbs)) {
          tempvect <- tempdbs[[cC]]
          tempvect <- tempvect[!is.na(tempvect)]
          if (length(tempvect)) {
            val <- if (bootstrap$average == "median") median(tempvect, na.rm=TRUE) else mean(tempvect, na.rm=TRUE)
            newstats[cR, (cC) := val]
          }
        }
      }
      
      # significance
      if (all(c("p.value","significance") %in% names(newstats))) {
        outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
        newstats[cR, significance := outPvalue$significance]
      }
      
      # bootstrap intervals
      if (all(c('F.value') %in% colnames(newstats))) {
        tempvect <- tempdbs[["F.value"]]
        tempvect <- tempvect[!is.na(tempvect)]
        ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
        if (!is.na(ci_boot[1])) {
          if (ci_boot[1] < 0) {ci_boot[1] <- 0}
          if (ci_boot[1] > newstats$F.value[cR]) {ci_boot[1] < newstats$F.value[cR]} 
        }
        if (!is.na(ci_boot[1])) {
          if (ci_boot[2] < 0) {ci_boot[2] <- 0}
          if (ci_boot[2] < newstats$F.value[cR]) {ci_boot[2] < newstats$F.value[cR]} 
        }
        newstats$F.value.ci.lower[cR] <- ci_boot[1]
        newstats$F.value.ci.upper[cR] <- ci_boot[2]
      }
      
      if (all(c('p.value') %in% colnames(newstats))) {
        tempvect <- tempdbs[["p.value"]]
        tempvect <- tempvect[!is.na(tempvect)]
        ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
        if (!is.na(ci_boot[1])) {
          if (ci_boot[1] < 0) {ci_boot[1] <- 0}
          if (ci_boot[1] > newstats$p.value[cR]) {ci_boot[1] < newstats$p.value[cR]} 
        }
        if (!is.na(ci_boot[1])) {
          if (ci_boot[2] < 0) {ci_boot[2] <- 0}
          if (ci_boot[2] < newstats$p.value[cR]) {ci_boot[2] < newstats$p.value[cR]} 
        }
        newstats$p.value.ci.lower[cR] <- ci_boot[1]
        newstats$p.value.ci.upper[cR] <- ci_boot[2]
      }
      
      if (all(c('ci.lower', 'ci.upper') %in% colnames(newstats))) {
        tempvect <- tempdbs[["effects"]]
        tempvect <- tempvect[!is.na(tempvect)]
        ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
        if (!is.na(ci_boot[1])) {
          if (ci_boot[1] < 0) {ci_boot[1] <- 0}
          if (ci_boot[1] > newstats$effects[cR]) {ci_boot[1] < newstats$effects[cR]} 
        }
        if (!is.na(ci_boot[1])) {
          if (ci_boot[2] < 0) {ci_boot[2] <- 0}
          if (ci_boot[2] < newstats$effects[cR]) {ci_boot[2] < newstats$effects[cR]} 
        }
        newstats$ci.lower[cR] <- ci_boot[1]
        newstats$ci.upper[cR] <- ci_boot[2]
      }
    }
    newstats <- as.data.frame(newstats)
    newstats$bootstrap <- TRUE
    
    # put data back
    textcall <- sprintf('results$%s <- as.data.frame(newstats)', statsofinterst[cSi])
    eval(parse(text=textcall))
    rm(newstats)
  }
  
  # remove folder
  unlink(bootstrap$tmpdir, recursive = TRUE)
  
  results <- lmerEffects2text(results)
  
  return(results)
}
