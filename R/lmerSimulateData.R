#' lmerSimulateData
#'
#' @description Function to create a simulated dataset from an lmer model. When method is conditionaldistribution, the function simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. The data is simulated multiple times and then subjects are removed until the desired sample sizes for each between subjects variable is obtained. When method is covariance, the function simulates data from a multivariate distribution using the covariance matrix using the MASS::mvrnorm function (when parametric is true) or using the mnonr::unonr function (when parametric is false).
#'
#' @param fit modelfit from the lmer function
#' @param between text list specifying the between subjects variables
#' @param within text list specifying any variables that should be considered as within subjects. If null assumes all factors are between subjects.
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label. If left NULL assumes all factors are between subjects
#' @param subsample decimal parameter 0 to 1 indicating the percent of the sample to use. default is the full sample. Ignored when method is conditionaldistribution
#' @param inflation decimal parameter indicating the group level multiplier for the final sample. default is no inflation
#' @param parametric boolean parameter indicating if the simulation should create a parametrically simulated dataset using the MASS::mvrnorm function (true) or a nonnormal simulated dataset using the mnonr::unonr function (false). Ignored when method is conditionaldistribution
#' @param method text specifying either conditionaldistribution (default) to use the simulate function or covariance to use the MASS::mvrnorm function or mnonr::unonr functions.
#'
#' @return results list output from the lmerEffects function
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 10, 2025
#'
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     newdata <- Rmimic::lmerSimulateData(fit, between=c('Group'), within=c('Time'), dependentvariable = 'Alertness', subjectid = 'PartID', inflation=1.0, method='conditionaldistribution')
#'     newdata2 <- Rmimic::lmerSimulateData(fit, between=c('Group'), within=c('Time'), dependentvariable = 'Alertness', subjectid = 'PartID', inflation=1.0, method='covariance', subsample=0.8)                
#'                
#' @importFrom dplyr slice_sample
#' @importFrom MASS mvrnorm
#' @importFrom stringr str_split
#' @importFrom mnonr mardia unonr
#' @importFrom stats model.frame formula cov runif
#' @importFrom data.table as.data.table data.table
#'
#' @export

lmerSimulateData <- function(fit, between=NULL, within=NULL, dependentvariable=NULL, subjectid=NULL, subsample=NULL, inflation=NULL, parametric=TRUE, method=NULL, screen=c(0.10, 0.90), constrain=NULL, verbose=FALSE) {
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(data.table)))))
  # Each repetition of the function would create a new simulated dataset
  # this approach should downweight any individual subject
  # simulation also allows for growing the sample
  
  if (!is.null(method)) {
    if (toupper(method) == toupper("covariance")) {
      method = "covariance"
    } else if (toupper(method) == toupper("parametric")) {
      method = "covariance"
      parametric = TRUE
    } else if (toupper(method) == toupper("nonparametric")) {
      method = "covariance"
      parametric = FALSE
    } else if (toupper(method) == toupper("conditionaldistribution")) {
      method = "conditionaldistribution"
    } else {
      method = "covariance"
    }
  } else {
    method = "covariance"
  }
  
  if (!is.null(constrain)) {
    if (toupper(constrain) == toupper("FALSE")) {
      constrain = FALSE
    } else if (toupper(constrain) == toupper("both")) {
      constrain = 'both'
    } else if (toupper(constrain) == toupper("upper")) {
      constrain = 'upper'
    } else if (toupper(constrain) == toupper("lower")) {
      constrain = 'lower'
    }
  } else {
    constrain = 'both'
  }
  
  if (method == "covariance") {
    # the dataset would be simulated based upon the covariance matrix of a subsample of data
    # covariance matrix is unique for each between subject level
    # random effects not including the participant id are considered as between subjects levels
    # covariance matrix is common for all within subject levels reflecting the inherent relationship of these levels
    
    if (verbose) {cat(sprintf('lmerSimulateData(): method: covariance\n'))}
    
    
    tempdbs <- data.table::as.data.table(stats::model.frame(fit))
    subject_ids <- unique(tempdbs[[subjectid]])
    randomfactor <- colnames(tempdbs)  
    randomfactor <- randomfactor[!(randomfactor %in% dependentvariable)]
    if (!is.null(within)) {
      within <- within[within %in% colnames(tempdbs)]
      randomfactor <- randomfactor[!(randomfactor %in% within)]
    }
    if (!is.null(between)) {
      between <- between[between %in% colnames(tempdbs)]
      randomfactor <- randomfactor[!(randomfactor %in% between)]
    }
    
    # Define number of subjects and observations per subject
    n_subjects <- length(subject_ids)
    
    # for each between subjects factor
    checknames <- colnames(tempdbs) # treat everything as between
    checknames <- checknames[which(checknames != dependentvariable)] # remove DV
    checknames <- checknames[which(checknames != subjectid)] # remove subjectID
    if (!is.null(within)) {
      #tempdbs$newvarnameforWT <- do.call(paste, c(tempdbs[within], sep="_by_"))
      tempdbs[, newvarnameforWT := do.call(paste, c(.SD, sep = "_by_")), .SDcols = within]
      checknames <- checknames[which(!(checknames %in% within))] # remove within
    } else {
      tempdbs$newvarnameforWT <- 'one'
    }
    if (length(between) > 0) {
      #tempdbs$newvarnameforBTW <- do.call(paste, c(tempdbs[between], sep="_by_"))
      tempdbs[, newvarnameforBTW := do.call(paste, c(.SD, sep = "_by_")), .SDcols = between]
      checknames <- checknames[which(!(checknames %in% between))] # remove between
    } else {
      tempdbs$newvarnameforBTW <- 'one'
    }
    uniqueBTW <- unique(tempdbs$newvarnameforBTW) 
    numberofsamples <- floor(n_subjects / length(uniqueBTW))
    if (!is.null(inflation)) {
      if (!is.numeric(inflation)) {
        inflation <- 1.0
      }
    }
    
    dataminval <- tempdbs[, min(get(dependentvariable), na.rm = TRUE)]
    datamaxval <- tempdbs[, max(get(dependentvariable), na.rm = TRUE)]
    
    mastertempdbs <- tempdbs
    
    # apply random effects
    randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
    randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
    randomformula <- randomformula[which(randomformula != dependentvariable)]
    randomformula <- stringr::str_split(randomformula, "\\+|\\*|\\/")[[1]]
    randomformula <- randomformula[which(randomformula != '')]
    
    if (verbose) {cat(sprintf('lmerSimulateData(): looping through each between subjects variable\n'))}
    
    BTWsimlist <- vector("list", length(uniqueBTW)) # Pre-allocate a list to store each iteration's result
    for (cBTW in 1:length(uniqueBTW)) {
      
      checkdata <- tempdbs[newvarnameforBTW == uniqueBTW[cBTW]]
      uniqueWT <- unique(checkdata$newvarnameforWT)
      uniqueIDs <- unique(checkdata[[subjectid]])
      numberofsamples <- floor(length(uniqueIDs))
      if (!is.null(inflation)) {
        numberofsamples <- floor(numberofsamples * inflation)
      }
      
      if (verbose) {cat(sprintf('lmerSimulateData(): creating covariance matrix for between subjects factor\n'))}
      
      # create table for correlation matrix - aggregate by subjectid and newvarnameforWT
      agg <- checkdata[, .(
                         meanDV = mean(get(dependentvariable), na.rm = TRUE),
                         sdDV   = sd(get(dependentvariable), na.rm = TRUE),
                         nDV    = .N
                       ), by = .(subject = get(subjectid), newvar = newvarnameforWT)
      ]
      outmatrix    <- data.table::dcast(agg, subject ~ newvar, value.var = "meanDV")
      outvarmatrix <- data.table::dcast(agg, subject ~ newvar, value.var = "sdDV")
      nmat         <- data.table::dcast(agg, subject ~ newvar, value.var = "nDV")
      outmatrix    <- as.matrix(outmatrix[, -1, with=FALSE])
      outvarmatrix <- as.matrix(outvarmatrix[, -1, with=FALSE])
      nmat         <- as.matrix(nmat[, -1, with=FALSE])
      samplecount <- rowMeans(nmat, na.rm = TRUE) # Compute samplecount (mean number of observations per subject)
      
      # transform to variance
      outvarmatrix <- outvarmatrix^2
      limits <- apply(outvarmatrix, 2, quantile, probs = screen, na.rm = TRUE)
      lower <- matrix(rep(limits[1, ], each = nrow(outvarmatrix)), nrow = nrow(outvarmatrix))
      upper <- matrix(rep(limits[2, ], each = nrow(outvarmatrix)), nrow = nrow(outvarmatrix))
      outvarmatrix[outvarmatrix < lower | outvarmatrix > upper] <- NA
      
      if (!is.null(subsample)) {
        if (is.numeric(subsample) && (subsample > 0) && (subsample < 1)) {
          subsample_fraction <- subsample
        } else {
          subsample_fraction <- 1  # No subsample
        }
        n <- nrow(outmatrix)
        if (n > 10) {
          subnumberofsamples <- floor(n * subsample_fraction)
        } else {
          subnumberofsamples <- n - 1
        }
        # Randomly select rows to keep
        if (subnumberofsamples < n) {
          keep_idx <- sample(n, subnumberofsamples, replace = FALSE)
          outmatrix    <- outmatrix[keep_idx, , drop = FALSE]
          outvarmatrix <- outvarmatrix[keep_idx, , drop = FALSE]
        }
      }
      
      if (length(uniqueWT) > 1) {
        tempsds <- colMeans(outvarmatrix, na.rm=TRUE)
      } else {
        tempsds <- mean(outvarmatrix, na.rm=TRUE)
      }
      
      # grab the means and prevent an issue
      if (length(uniqueWT) > 1) {
        tempmeans <- colMeans(outmatrix, na.rm=TRUE)
        tempsigma <- tryCatch(
          {stats::cov(outmatrix, use = "pairwise.complete.obs")},
          error = function(e) {
            # stupid simple hack for insufficient observations
            outmatrix <- rbind(rbind(outmatrix, outmatrix), outmatrix)
            tempsigma <- stats::cov(outmatrix, use="pairwise.complete.obs")
          }
        )
      } else {
        tempmeans <- mean(outmatrix, na.rm=TRUE)
        tempsigma <- matrix(c(tempsds,1),1,1)
      }
      
      multiplecasespersubject <- mean(samplecount, na.rm=TRUE, trim=0.25)
      simlist <- vector("list", multiplecasespersubject) # Pre-allocate a list to store each iteration's result
      for (cMCpS in 1:multiplecasespersubject) {
        if (verbose) {cat(sprintf('lmerSimulateData(): simulating multivariate data for case %d\n', cMCpS))}
        if (length(uniqueWT) > 1) {
          
          # Simulate Multivariate Data
          if (parametric) {
            withindata <- invisible(suppressMessages(suppressWarnings(MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma))))
          } else {
            skewkurtcheck <- invisible(mnonr::mardia(outmatrix))
            skewkurtcheckuni <- data.table::as.data.table(skewkurtcheck$univariate)
            rownames(skewkurtcheckuni) <- rownames(skewkurtcheck$univariate)
            withindata <- invisible(suppressWarnings(mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=unlist(skewkurtcheckuni[1,]), kurtosis = unlist(skewkurtcheckuni[3,]))))
            colnames(withindata) <- colnames(outmatrix)
          }
          withindata <- data.table::as.data.table(withindata)
          
          # unpack data
          dt_long <- data.table::data.table(
            newvarnameforWT = rep(uniqueWT, each = numberofsamples),
            newvarnameforBTW = rep(uniqueBTW[cBTW], numberofsamples * length(uniqueWT)),
            subjectid = sprintf('IDCode_%s_%d', uniqueBTW[cBTW], rep(seq_len(numberofsamples), times = length(uniqueWT)))
          )
          data.table::setnames(dt_long, "subjectid", subjectid)
          dt_long[, (dependentvariable) := unlist(withindata, use.names = FALSE)]
          
          # Add missing columns as NA to match tempdbs
          missing_cols <- setdiff(colnames(tempdbs), colnames(dt_long))
          if (length(missing_cols) > 0) {
            dt_long[, (missing_cols) := NA]
          }
          data.table::setcolorder(dt_long, colnames(tempdbs))
          simlist[[cMCpS]] <- dt_long  # Store in list
          
        } else {
          # Simulate Univariate Data
          if (parametric) {
            withindata <- invisible(suppressMessages(suppressWarnings(MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma))))
          } else {
            withindata <- invisible(suppressMessages(suppressWarnings(mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=NULL, kurtosis = NULL))))
            colnames(withindata) <- colnames(outmatrix)
          }
          withindata <- data.table::as.data.table(withindata)
          
          # unpack data
          dt_long <- data.table::data.table(
            newvarnameforWT = rep(uniqueWT, each = numberofsamples),
            newvarnameforBTW = rep(uniqueBTW[cBTW], numberofsamples * length(uniqueWT)),
            subjectid = sprintf('IDCode_%s_%d', uniqueBTW[cBTW], rep(seq_len(numberofsamples), times = length(uniqueWT)))
          )
          data.table::setnames(dt_long, "subjectid", subjectid)
          dt_long[, (dependentvariable) := unlist(withindata, use.names = FALSE)]
          
          # Add missing columns as NA to match tempdbs
          missing_cols <- setdiff(colnames(tempdbs), colnames(dt_long))
          if (length(missing_cols) > 0) {
            dt_long[, (missing_cols) := NA]
          }
          data.table::setcolorder(dt_long, colnames(tempdbs))
          simlist[[cMCpS]] <- dt_long  # Store in list
        }
      } # cases per subject
      BTWsimlist[[cBTW]] <- data.table::rbindlist(simlist, use.names = TRUE, fill = TRUE)
      
    } # between subjects
    mastersim_data <- data.table::rbindlist(BTWsimlist, use.names = TRUE, fill = TRUE)
    
    # decode
    if (length(between) > 0) {
      mastersim_data[, (between) := data.table::tstrsplit(newvarnameforBTW, '_by_')]
    }
    if (!is.null(within)) {
      mastersim_data[, (within) := data.table::tstrsplit(newvarnameforWT, '_by_')]
    }
    
    # constrain data to observed values
    if ((constrain == 'lower') | (constrain == 'both')) {
      mastersim_data[get(dependentvariable) < dataminval, (dependentvariable) := dataminval]
    }
    if ((constrain == 'upper') | (constrain == 'both')) {
      mastersim_data[get(dependentvariable) > datamaxval, (dependentvariable) := datamaxval]
    }
    
    # restore unadjusted original data - just in case
    tempdbs <- mastertempdbs
    
    # Apply original precision
    try({
      vals <- tempdbs[[dependentvariable]]
      vals <- vals[!is.na(vals)]
      get_decimals <- function(x) {ifelse(is.na(x), 0, nchar(sub("^\\d+\\.?|0+$", "", formatC(x, digits=20, format="f")))) }
      mastersim_data[, (dependentvariable) := round(get(dependentvariable), digits = max(get_decimals(vals)))]
    }, silent = TRUE)
    
    
    
    # populate other variables
    remainingfactors <- randomfactor
    remainingfactors <- remainingfactors[which(remainingfactors != subjectid)] # remove subjectID
    if (length(remainingfactors) > 0) {
      for (cRF in 1:length(remainingfactors)) {
        # does this factor ever vary within a participant
        withinfactorRF <- tempdbs[, .(n_unique = data.table::uniqueN(get(remainingfactors[cRF]))), by = subjectid][, any(n_unique > 1)]
        
        if (!(withinfactorRF)) {
          # between subjects factor
          
          subworkingdatabase <- tempdbs[, .(N = .N), by = c(subjectid, remainingfactors[cRF])]
          uniqueRF <- unique(subworkingdatabase[[remainingfactors[cRF]]])
          
          # allocation ratio per group level
          alloc_counts <- subworkingdatabase[, .(N = .N), by = eval(remainingfactors[cRF])]
          alloc_counts[, ratio := N / sum(N)]
          alloc_counts[, N := NULL]  # remove count, keep ratio
          totaloutsubjects <- data.table::uniqueN(mastersim_data[[subjectid]])
          alloc_counts[, n_subjects := round(ratio * totaloutsubjects)]
          diff <- totaloutsubjects - sum(alloc_counts$n_subjects)
          if (diff != 0) alloc_counts$n_subjects[1] <- alloc_counts$n_subjects[1] + diff
          allocdistr <- unlist(mapply(rep, alloc_counts[[remainingfactors[cRF]]], alloc_counts$n_subjects))
          allocdistr <- sample(allocdistr)  # shuffle
          uniqueoutsubjects <- unique(mastersim_data[[subjectid]])
          
          dt_alloc <- data.table::data.table(subid = uniqueoutsubjects, allocdistr = allocdistr)
          data.table::setnames(dt_alloc, "subid", subjectid)
          data.table::setnames(dt_alloc, "allocdistr", remainingfactors[cRF])
          idx <- match(mastersim_data[[subjectid]], dt_alloc[[subjectid]])
          mastersim_data[, (remainingfactors[cRF]) := dt_alloc[[remainingfactors[cRF]]][idx]]
        } # is a between subjects factor
      } # each remaining factor
    }
    
     # remove decode
     mastersim_data$newvarnameforBTW <- NULL # remove it
     mastersim_data$newvarnameforWT <- NULL # remove it 
     
     # swap ID code for simpler form
     unique_ids <- unique(mastersim_data[[subjectid]])
     replacement_codes <- sprintf("ID%02d", seq_along(unique_ids))
     map <- setNames(replacement_codes, unique_ids)
     mastersim_data[, (subjectid) := map[ mastersim_data[[subjectid]] ] ]
     
    # return same order
    tempdbs <- stats::model.frame(fit)
    data.table::setcolorder(mastersim_data, colnames(tempdbs))
    mastersim_data <- mastersim_data[, .SD, .SDcols = colnames(tempdbs)]
    mastersim_data <- as.data.frame(mastersim_data)
    
  } else {
    #method = "conditionaldistribution"
    
    if (is.null(inflation)) {
      inflation <- 1.0
    } else {
      if (!is.numeric(inflation)) {
        inflation <- 1.0
      }
    }
    
    tempdbs <- stats::model.frame(fit)
    if (length(between) > 0) {
      tempdbs$newvarnameforBTW <- do.call(paste, c(tempdbs[between], sep="_by_"))
    } else {
      tempdbs$newvarnameforBTW <- 'one'
    }
    uniqueids <- unique(tempdbs[,subjectid])
    uniqueBTW <- unique(tempdbs$newvarnameforBTW)
    uniqueBTWL <- rep_len(NA, length(uniqueBTW))
    for (cBTW in 1:length(uniqueBTW)) {
      currentcount <- 0
      for (cSID in 1:length(uniqueids)) {
        checkindx <- which(tempdbs[,subjectid] == uniqueids[cSID] & tempdbs$newvarnameforBTW == uniqueBTW[cBTW])
        if (length(checkindx) > 0) {
          currentcount <- currentcount + 1
        }
      }
      uniqueBTWL[cBTW] <- currentcount
    }
    numberofsamples <- floor(inflation * uniqueBTWL)
    
    tempdbs <- data.table::as.data.table(tempdbs)
    dataminval <- tempdbs[, min(get(dependentvariable), na.rm = TRUE)]
    datamaxval <- tempdbs[, max(get(dependentvariable), na.rm = TRUE)]
    
    mastersim_data <- invisible(lmerSimulate_conditionaldistribution(fit, dependentvariable=dependentvariable, subjectid=subjectid, between=between, targetN=numberofsamples, subsample=subsample))
    mastersim_data <- data.table::as.data.table(mastersim_data)
    
    # swap ID code for simpler form
    unique_ids <- unique(mastersim_data[[subjectid]])
    replacement_codes <- sprintf("ID%02d", seq_along(unique_ids))
    map <- setNames(replacement_codes, unique_ids)
    mastersim_data[, (subjectid) := map[ mastersim_data[[subjectid]] ] ]
    
    # constrain data to observed values
    if ((constrain == 'lower') | (constrain == 'both')) {
      mastersim_data[get(dependentvariable) < dataminval, (dependentvariable) := dataminval]
    }
    if ((constrain == 'upper') | (constrain == 'both')) {
      mastersim_data[get(dependentvariable) > datamaxval, (dependentvariable) := datamaxval]
    }
    mastersim_data <- as.data.frame(mastersim_data)
    
  }
  return(mastersim_data)
}





  
