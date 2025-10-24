#' lmerEffectsBootstrapContrast
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerPosthoc function. This portion handles the test contrasts.
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

lmerEffectsBootstrapContrast <- function(results, contrastin, bootstrap, ...) {
  
  dots <- tryCatch({
    dots <- list(...) 
  }, error = function(e) {
    dots <- list() 
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
  
  # Define function for MAD outliers
  runmad <- function(x) {
    idx_outliers <- integer(0)
    # ensure numeric
    if (!is.numeric(x)) {
      return(x)
    }
    x_ok <- x[!is.na(x) & is.finite(x)]
    if (length(x_ok) >= 3) {
      med <- tryCatch(
        stats::median(x_ok),
        error = function(e) {
          #warning("Median computation failed: ", conditionMessage(e))
          return(NA_real_)
        }
      )
      mad_raw <- tryCatch(
        stats::mad(x_ok, constant = 1, na.rm = TRUE),
        error = function(e) {
          #warning("MAD computation failed: ", conditionMessage(e))
          return(NA_real_)
        }
      )
      # scale MAD to be comparable to SD (consistent with normal): *1.4826
      scale <- tryCatch({
        if (is.finite(mad_raw) && mad_raw > 0) {
          mad_raw * 1.4826
        } else {
          stats::IQR(x_ok, na.rm = TRUE) / 1.349
        }
      }, error = function(e) {
        #warning("Scale computation failed: ", conditionMessage(e))
        return(NA_real_)
      })
      if (is.finite(scale) && scale > 0 && !is.na(med)) {
        z <- tryCatch(
          abs(x - med) / scale,
          error = function(e) {
            #warning("Z-score computation failed: ", conditionMessage(e))
            return(rep(NA_real_, length(x)))
          }
        )
        idx_outliers <- which(z > 3.5) # MAD threshold for outliers (robust)
      }
      # replace outliers with NA, then drop them
      if (length(idx_outliers) > 0) {
        x[idx_outliers] <- NA
        x <- x[!is.na(x)]
      }
    }
    return(x)
  }
  
  options(warn = -1)
  
  # Define function for a single repetition
  run_two <- function(results, resample_min, resample_max, subsample, inflation, method) {
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
      smp <- invisible(suppressWarnings(suppressMessages(lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=TRUE, method = "covariance"))))
    }
    if (method == "nonparametric") {
      # Useful for growing the sample 
      # keeps the exact same effect size by allowing the group variation to increase with larger samples
      smp <- invisible(suppressWarnings(suppressMessages(lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, parametric=FALSE, method = "covariance"))))
    }
    if (method == "default") {
      # works the best as it is a wrapper around simulate - not ideal for growing the sample as the effect size will grow with it but will keep the data around the original mean and standard deviation
      smp <- invisible(suppressWarnings(suppressMessages(lmerSimulateData(results$fit, between=c(results$between, results$covariates), within=results$within, dependentvariable=results$dependentvariable, subjectid=results$subjectid, subsample=subsample, inflation=inflation, method = "conditionaldistribution"))))
    }
    
    if (!is.null(smp)) {
      smp <- as.data.frame(smp)
      # rerun model on new data
      newfit <- tryCatch({
        #newfit <- invisible(suppressWarnings(suppressMessages(update(results$fit, data=smp, evaluate = TRUE))))
        newfit <- invisible(suppressWarnings(suppressMessages(update(results$fit, formula = formula(results$fit), data = smp, evaluate = TRUE))))
      }, error = function(e) {
        cat(sprintf('lmerSimulateData - model failure\n'))
        newfit <- NULL
      })
    } else {
      newfit <- NULL
    }
    
    if (!is.null(newfit)) {
      return(newfit)
    } else {
      return(NULL)
    }
  }
  
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future)))))
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(future.apply)))))
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(progressr)))))

  totalsample <- nrow(stats::model.frame(results$fit))
  if (is.null(bootstrap$resample_min)) {
    bootstrap$resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (is.null(bootstrap$resample_max)) {
    bootstrap$resample_max <- nrow(stats::model.frame(results$fit))
  }  
  
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
    #tmpdir <- file.path(getwd(), paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), collapse = ""))
    bootstrap$tmpdir <- tempfile(pattern = "anovacontrast", tmpdir = getwd(), fileext = "")
  }
  dir.create(bootstrap$tmpdir, showWarnings = FALSE)
  
  # Check how many result files exist
  file_list <- list.files(bootstrap$tmpdir, pattern = "^result_contrast_.*\\.csv$")
  
  # Update to reflect if any files have already been created
  repetitions <- repetitions - length(file_list) + 1
  
  if (length(file_list) < max_files) {
    
    # Enable progress handler
    #handlers("txtprogressbar") 
    handlers(list(
      handler_progress(
        format   = "lmerEffectsBootstrapContrasts() simulating data [:bar] :percent :eta",
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
        file_list <- list.files(bootstrap$tmpdir, pattern = "^result_contrast_.*\\.csv$")
        
        if (length(file_list) < max_files) {
          smp <- tryCatch({
            simfit <- run_two(results, bootstrap$resample_min, bootstrap$resample_max, bootstrap$subsample, bootstrap$inflation, bootstrap$method)
            if (!is.null(simfit)) {
              outtable <- lmerPosthocsubprocessContrasts(simfit, emmeansdf=contrastin$emmeansdf, effectofinterest=contrastin$effectofinterest, currentfactor=contrastin$currentfactor,
                                                         factorsinvolved=contrastin$factorsinvolved, otherfactorsinvolved=contrastin$otherfactorsinvolved,
                                                         decomptext=contrastin$decomptext, factortag=contrastin$factortag,
                                                         dependentvariable=results$dependentvariable[1], 
                                                         subjectid=results$subjectid[1], confidenceinterval=results$confidenceinterval,
                                                         studywiseAlpha=results$studywiseAlpha, within=results$within)
              data.table::fwrite(outtable, tempfile(pattern = "result_contrast_", tmpdir = bootstrap$tmpdir, fileext = ".csv"), append = FALSE)
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
  file_list <- list.files(bootstrap$tmpdir, pattern = "^result_contrast_.*\\.csv$", full.names=TRUE)
  if (length(file_list) > 100) {
    #datain <- data.table::rbindlist(future_lapply(file_list, data.table::fread), use.names=TRUE, fill=TRUE)
    # Enable progress handler
    #handlers("txtprogressbar") 
    handlers(list(
      handler_progress(
        format   = "lmerEffectsBootstrapContrasts() compiling data [:bar] :percent :eta",
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
  newstats <- data.table::as.data.table(data.table::copy(contrastin$outtable))
  
  # --- summary stats ---
  cols1 <- c("df","t.ratio","p.value","effectsize","correlation")
  agg1 <- datain[, lapply(.SD, function(x) {
    x <- runmad(x) 
    if (bootstrap$average=="median") median(x, na.rm=TRUE) else mean(x, na.rm=TRUE)
  }), by=.(idtag, hold, contrast, decomp), .SDcols=cols1]
  
  # --- confidence intervals ---
  ci_fun <- function(x) {
    x <- x[!is.na(x)]
    x <- runmad(x) 
    if (length(x)==0) return(c(NA_real_, NA_real_))
    qs <- quantile(x, c(.05,.95), na.rm=TRUE)
    c(lower=max(0, qs[1], na.rm=TRUE), upper=qs[2])
  }
  vars <- c("effectsize","t.ratio","p.value")
  ci_wide <- datain[, {
    vals <- unlist(lapply(vars, function(v) ci_fun(get(v))))
    nm   <- unlist(lapply(vars, function(v) paste0(v, c(".conf.int.lower",".conf.int.upper"))))
    setNames(as.list(vals), nm)
  }, by = .(idtag, hold, contrast, decomp)]
  # column-bind all CI columns to agg1
  data.table::setkey(ci_wide, idtag, hold, contrast, decomp)
  agg1 <- ci_wide[agg1, on = .(idtag, hold, contrast, decomp)]
  
  data.table::setnames(agg1, old='t.ratio.conf.int.upper', new='t.conf.int.upper')
  data.table::setnames(agg1, old='t.ratio.conf.int.lower', new='t.conf.int.lower')
  data.table::setnames(agg1, old='p.value.conf.int.lower', new='p.conf.int.lower')
  data.table::setnames(agg1, old='p.value.conf.int.upper', new='p.conf.int.upper')
  
  cols2 <- c("C1n","C1mean","C1sd","C2n","C2mean","C2sd")
  agg2 <- datain[, lapply(.SD, function(x) {
    x <- runmad(x) 
    if (bootstrap$average=="median") median(x, na.rm=TRUE) else mean(x, na.rm=TRUE)
  }), by=.(idtag, hold, contrast, decomp), .SDcols=cols2]
  # column-bind all columns to agg1
  data.table::setkey(agg2, idtag, hold, contrast, decomp)
  agg1 <- agg2[agg1, on = .(idtag, hold, contrast, decomp)]
  agg1[, significant := NA]
  
  agg3 <- datain[,.(C1name <- unique(C1name),C2name <- unique(C2name)),
                 by=.(idtag, hold, contrast, decomp)]
  data.table::setnames(agg3, old='V1', new='C1name')
  data.table::setnames(agg3, old='V2', new='C2name')
  data.table::setkey(agg3, idtag, hold, contrast, decomp)
  agg1 <- agg3[agg1, on = .(idtag, hold, contrast, decomp)]
  
  data.table::setcolorder(agg1, names(newstats))
  
  agg1 <- as.data.frame(agg1)
  newstats <- as.data.frame(newstats)
  # pull significance
  for (cdtR in 1:nrow(agg1)) {
    outPvalue <- fuzzyP(agg1$p.value[cdtR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
    agg1$significant[cdtR] <- outPvalue$significance
  }
  
  if (bootstrap$reporteddata == "actual") {
    # make sure we didnt get a name position swap
    checkonce <- TRUE
    for (cR in 1:nrow(newstats)) {
      tempdbs <- agg1[which(agg1$idtag == newstats$idtag[cR] & agg1$hold == newstats$hold[cR]),]
      if (nrow(tempdbs) == 0) {
        # data simulation likely swapped position of c1 and c2
        tempvect <- stringr::str_split(newstats$idtag[cR], "-")[[1]]
        newtempvect <- tempvect
        newtempvect[length(tempvect)] <- tempvect[length(tempvect)-2]
        newtempvect[length(tempvect)-2] <- tempvect[length(tempvect)]
        newstats$idtag[cR] <- paste0(newtempvect, collapse="-")
        if (checkonce) {
          c1name <- newstats$C1name
          c2name <- newstats$C2name
          newstats$C1name <- c2name
          newstats$C2name <- c1name
          
          c1n <- newstats$C1n
          c2n <- newstats$C2n
          newstats$C1n <- c2n
          newstats$C2n <- c1n
          
          c1n <- newstats$C1mean
          c2n <- newstats$C2mean
          newstats$C1mean <- c2n
          newstats$C2mean <- c1n
          
          c1n <- newstats$C1sd
          c2n <- newstats$C2sd
          newstats$C1sd <- c2n
          newstats$C2sd <- c1n
          checkonce <- FALSE
        }
      }
    }
    # swap values with actual data
    agg1 <- data.table::as.data.table(agg1)
    newstats <- data.table::as.data.table(newstats)
    keycols <- c("idtag","hold","contrast","decomp")
    cols    <- c("C1n", "C2n","C1mean","C2mean","C1sd","C2sd")
    agg1[newstats, (cols) := mget(paste0("i.", cols)), on = .(idtag, hold, contrast, decomp)]
    agg1 <- as.data.frame(agg1)
    newstats <- as.data.frame(newstats)
  }
  agg1$bootstrap <- TRUE
  outtable <- agg1 
  
  # remove folder
  unlink(bootstrap$tmpdir, recursive = TRUE)
  
  return(outtable)
}
