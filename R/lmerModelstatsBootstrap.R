#' lmerModelstatsBootstrap
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerModelstats function. 
#'
#' @param results List containing output from lmerModelstats
#' @param repetitions Numeric entry indicating the number of repetitions to perform
#' @param between text list specifying the between subjects variables and any between subjects covariates
#' @param within text list specifying any variables that should be considered as within subjects. If null assumes all variables are between subjects.
#' @param subjectid text specifying the subject id label.
#' @param method String entry indicating the approach. Default simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. Resample performs data resampling with replacement from the existing dataset. Parametric simulates data from a multivariate normal distribution using the MASS::mvrnorm function. Nonparametric simulates data from a multivariate nonnormal distribution using the mnonr::unonr function.
#' @param subsample Numeric entry 0 to 1 indicating what percentage of the data to use for informing the parametric simulation. Ignored for method resample.
#' @param inflation Numeric entry indicating how many times the original sample to simulate too. Ignored for method resample.
#' @param resample_min Numeric entry indicating the minimum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param resample_max Numeric entry indicating the maximum number of samples to include (with replacement). Default keeps the total samples the same. Ignored for method parametric.
#' @param average text parameter to indicate how the results should be collapsed. Default is using the median. Other options is mean.
#' @param progressbar Boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.

#' @return
#' \item{results}{List containing boostrapped output.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, June 4, 2025
#'
#' @importFrom dplyr slice_sample
#' @importFrom progress progress_bar
#' @importFrom Rmisc CI
#' @importFrom miscTools colMedians
#' @importFrom stats runif model.frame
#' 
#' @examples
#' \dontrun{
#'     altfit <- lmerTest::lmer(Alertness ~ Time + (1 | PartID), data = Rmimic::alertness)
#'     fit <- lmerTest::lmer(Alertness ~ Group + Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerModelstats(fit, altfit)
#'     results <- Rmimic::lmerModelstatsBootstrap(results, repetitions=999, between=c('Group'),
#'                                               within=c('Time'), subjectid = "PartID",
#'                                               subsample=0.96, method='default', progressbar=TRUE)
#'     }
#'
#' @export

lmerModelstatsBootstrap <- function(results, repetitions, between, within=NULL, subjectid=NULL, resample_min=NULL, resample_max=NULL, subsample=0.96, inflation=1.0, method='default', average='median', progressbar=TRUE) {
  
  
  totalsample <- nrow(stats::model.frame(results$fit))
  if (is.null(resample_min)) {
    resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (is.null(resample_max)) {
    resample_max <- nrow(stats::model.frame(results$fit))
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
  
  # CHECK TO MAKE SURE THAT lmerEffects HAS BEEN RUN
  listofoutputsfromlmerEffects <- c('fit', 'Chisq', 'deviance', 'coefficients')
  if (all(listofoutputsfromlmerEffects %in% names(results))) {
    
    if (progressbar) {
      # establish progress
      #cat(sprintf('lmerModelstatsBootstrap() beginning processing '))
      countticks <- 0
      stepcounts <- 10
      checkins <- floor(seq(0, repetitions, length.out = stepcounts+1))
      pb <- progress::progress_bar$new(total = stepcounts,
                                       format = " lmerModelstatsBootstrap() processing [:bar] :percent eta: :eta",
                                       clear = TRUE, width= 120)
      pb$tick(0)
    }
    
    workingres <- NULL
    
    resstore <- list()
    
    # loop through analysis
    for (cN in 1:repetitions) {
      # populates dataset subsampled from the original sample with replacement
      if (method == "resample") {
        # Determine how large a random sample
        numberofsamples <- floor(stats::runif(1, min=resample_min, max=resample_max))
        smp <- dplyr::slice_sample(stats::model.frame(results$fit), n=numberofsamples, replace=TRUE)
      }
      if (method == "parametric") {
        # Useful for growing the sample 
        smp <- Rmimic::lmerSimulateData(results$fit, between=between, within=within, dependentvariable=results$dv, subjectid=subjectid, subsample=subsample, inflation=inflation, parametric=TRUE, method = "covariance")
      }
      if (method == "nonparametric") {
        # Useful for growing the sample 
        smp <- Rmimic::lmerSimulateData(results$fit, between=between, within=within, dependentvariable=results$dv, subjectid=subjectid, subsample=subsample, inflation=inflation, parametric=FALSE, method = "covariance")
      }
      if (method == "default") {
        # works the best as it is a wrapper around simulate - not ideal for growing the sample as the effect size will grow with it
        smp <- Rmimic::lmerSimulateData(results$fit, between=between, within=within, dependentvariable=results$dv, subjectid=subjectid, subsample=subsample, inflation=inflation, method = "conditionaldistribution")
      }
      
      # rerun model on new data
      newfit <- tryCatch({
        newfit <- update(results$fit, data=smp, evaluate = TRUE)
      }, error = function(e) {
        cat(sprintf('lmerModelstatsBootstrap - model failure\n'))
        newfit <- NULL
      })
      
      if (!is.null(newfit)) {
        # update the alternative model if needed
        usablerun <- TRUE
        if (!is.null(results$altfit)) {
          newaltfit <- tryCatch({
            newaltfit <- update(results$altfit, data=smp, evaluate = TRUE)
          }, error = function(e) {
            cat(sprintf('lmerModelstatsBootstrap - model failure\n'))
            usablerun <- FALSE
            newaltfit <- NULL
          })
          if (is.null(newaltfit)) {
            usablerun <- FALSE
          }
        }
        if (usablerun) {
          # extract information
          newresults <- lmerModelstats(newfit, newaltfit, df = results$dfparam, confidenceinterval=results$confidenceinterval, studywiseAlpha=results$studywiseAlpha)
          newresults$numparticipants <- length(unique(smp[,subjectid]))
          
          # store it
          textcall <- sprintf('resstore$repetition%d <- newresults', cN)
          eval(parse(text=textcall))
        }
        
      }
      if (progressbar) {
        # update progress
        if (cN %in% checkins) {
          #  cat(sprintf('.'))
          pb$tick()
          countticks <- countticks + 1
        }
      }
    }
    if (progressbar) {
      #cat(sprintf('\n'))
      if (countticks < stepcounts) {
        pb$tick(stepcounts)
      }
    }
    
    # summarize
    results <- resmerge(results, resstore, average=average)
    
    # obtain text outputs
    results <- lmerModelstats2text(results, confidenceinterval=results$confidenceinterval, studywiseAlpha=results$studywiseAlpha, testconfidence=TRUE, significanceconfidence=TRUE)
    
    
    # tag messageout 
    if ('messageout' %in% names(results)) {
      
      outstring <- sprintf('Unstandardized effects were computed based upon bootstrapped analyses with')
      if (method == "resample") {
        outstring <- sprintf('%s %d resamples', outstring, repetitions)
        if (!((resample_min == nrow(stats::model.frame(results$fit))) & (resample_max == nrow(stats::model.frame(results$fit))))) {
          outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%)',
                               outstring, resample_min, round((resample_min/totalsample)*100, digits=1),
                               resample_max, round((resample_max/totalsample)*100, digits=1))
        }
        outstring <- sprintf('%s of the original dataset (with replacement).', outstring)
      }
      if ((method == "parametric") | (method == "nonparametric")) {
        if (method == "parametric") {
          outstring <- sprintf('%s %d datasets simulated from a multivariate normal distribution using the MASS mvrnorm function (Venables &#38; Ripley, 2002).', outstring, repetitions)
        } else {
          outstring <- sprintf('%s %d datasets simulated from a multivariate non-normal distribution using the mnonr unonr function (Qu &#38; Zhang, 2020).', outstring, repetitions)
        }
        outstring <- sprintf('%s For each simulation the covariance matrix was informed by', outstring)
        if (subsample < 1.0) {
          outstring <- sprintf('%s a subsample of %.1f%% of the original data', outstring, round((subsample)*100, digits=1))
        } else {
          outstring <- sprintf('%s the full sample of original data', outstring)
        }
        if (inflation > 1.0) {
          outstring <- sprintf('%s and extrapolated to a final sample of %.0f participants', outstring, round(results$meanofparticipants, digits=0))
        } else {
          outstring <- sprintf('%s.', outstring)
        }
      }
      if (method == "default") {
        outstring <- sprintf('%s %d datasets simulated by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters', outstring, repetitions)
        
        if (!((resample_min == nrow(stats::model.frame(results$fit))) & (resample_max == nrow(stats::model.frame(results$fit))))) {
          outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%) of the original dataset (with replacement)',
                               outstring, resample_min, round((resample_min/totalsample)*100, digits=1),
                               resample_max, round((resample_max/totalsample)*100, digits=1))
        }
        outstring <- sprintf('%s.', outstring)
        
        if (!is.null(subsample)) {
          if (subsample < 1.0) {
            outstring <- sprintf('%s For each simulation the the conditional distribution of the outcome variable was informed by', outstring)
            outstring <- sprintf('%s a subsample of %.1f%% of the original data within each between subjects factor.', outstring, round((subsample)*100, digits=1))
          }
        }
        
      }
      
      results$messageout <- sprintf('%s %s', results$messageout, outstring)
    }
    
  }
  
  return(results)
}


resmerge <- function(results, resstore, average='median') {
  
  varsofinterest <- c("AIC", "df", "Chisq", "F.value", "deviance", "p.value", "VIF", "fixed.r.squared", "random.r.squared", "model.r.squared","fsquared","fsquared.ci.lower", "fsquared.ci.upper")
  for (cR in 1:length(varsofinterest)) {
    newresults <- unlist(resresearch(resstore, target=varsofinterest[cR])) # should return a list
    if (average == 'median') {
      newresults <- median(newresults, na.rm=TRUE)
    } else {
      newresults <- mean(newresults, na.rm=TRUE)
    }
    textcall <- sprintf('results$%s <- newresults[1]', varsofinterest[cR])
    eval(parse(text=textcall))
  }
  
  # significance test
  outPvalue <- Rmimic::fuzzyP(results$p.value, studywiseAlpha=results$studywiseAlpha, html=TRUE)
  results$significant <- outPvalue$significance
  
  # additional confidence intervals
  newresults <- unlist(resresearch(resstore, target='Chisq')) # should return a list
  civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
  results$Chisq.ci.lower <- civals[which(names(civals) == 'lower')]
  results$Chisq.ci.upper <- civals[which(names(civals) == 'upper')]
  
  newresults <- unlist(resresearch(resstore, target='p.value')) # should return a list
  civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
  results$p.value.ci.lower <- civals[which(names(civals) == 'lower')]
  if (results$p.value.ci.lower < 0) {results$p.value.ci.lower <- 0}
  results$p.value.ci.upper <- civals[which(names(civals) == 'upper')]
  if (results$p.value.ci.upper < 0) {results$p.value.ci.upper <- 0}
  
  newresults <- unlist(resresearch(resstore, target='fixed.r.squared')) # should return a list
  civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
  results$fixed.r.squared.ci.lower <- civals[which(names(civals) == 'lower')]
  results$fixed.r.squared.ci.upper <- civals[which(names(civals) == 'upper')]
  
  newresults <- unlist(resresearch(resstore, target='random.r.squared')) # should return a list
  civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
  results$random.r.squared.ci.lower <- civals[which(names(civals) == 'lower')]
  results$random.r.squared.ci.upper <- civals[which(names(civals) == 'upper')]
  
  
  
  # Coefficients
  newresults <- resresearch(resstore, target='coefficients') # should return a table
  results$coefficients$significant <- NA
  results$coefficients$text <- NA
  # create empty variables
  results$coefficients$p.value.ci.lower <- NA
  results$coefficients$p.value.ci.upper <- NA
  results$coefficients$t.ci.lower <- NA
  results$coefficients$t.ci.upper <- NA
  for (cR in 1:nrow(results$coefficients)) {
    tempdbs <- newresults[which(newresults$Variable == results$coefficients$Variable[cR]),]
    colsofinterest <- c("B", "B.lower.conf.int", "B.upper.conf.int", "SE", "Beta", "Beta.lower.conf.int", "Beta.upper.conf.int",
                        "t", "df", "p.value")
    for (cC in 1:length(colsofinterest)) {
      if (average == 'median') {
        results$coefficients[cR,colsofinterest[cC]] <- median(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      } else {
        results$coefficients[cR,colsofinterest[cC]] <- mean(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      }
    }
    # significance test
    outPvalue <- Rmimic::fuzzyP(results$coefficients$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
    results$coefficients$significant[cR] <- outPvalue$significance
    
    civals <- suppressWarnings(Rmisc::CI(tempdbs$p.value, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$coefficients$p.value.ci.lower[cR] <- civals[which(names(civals) == 'lower')]
    if (results$coefficients$p.value.ci.lower[cR] < 0) {results$coefficients$p.value.ci.lower[cR] <- 0}
    results$coefficients$p.value.ci.upper[cR] <- civals[which(names(civals) == 'upper')]
    if (results$coefficients$p.value.ci.upper[cR] < 0) {results$coefficients$p.value.ci.upper[cR] <- 0}
    
    
    civals <- suppressWarnings(Rmisc::CI(tempdbs$t, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$coefficients$t.ci.lower[cR] <- civals[which(names(civals) == 'lower')]
    results$coefficients$t.ci.upper[cR] <- civals[which(names(civals) == 'upper')]
    
  }
  
  # change statistics
  if (!is.null(results$altfit)) {
    
    results$change$significant <- NA
    results$changetext <- NA
    varsofinterest <- c("AIC", "df", "Chisq", "F.value", "deviance", "p.value", "VIF", "fixed.r.squared", "random.r.squared", "model.r.squared","fsquared","fsquared.ci.lower", "fsquared.ci.upper")
    for (cR in 1:length(varsofinterest)) {
      newresults <- unlist(resresearch(resstore, target=sprintf('change$%s', varsofinterest[cR]))) # should return a list
      if (average == 'median') {
        newresults <- median(newresults, na.rm=TRUE)
      } else {
        newresults <- mean(newresults, na.rm=TRUE)
      }
      textcall <- sprintf('results$change$%s <- newresults[1]', varsofinterest[cR])
      eval(parse(text=textcall))
    }
    
    # significance test
    outPvalue <- Rmimic::fuzzyP(results$change$p.value, studywiseAlpha=results$studywiseAlpha, html=TRUE)
    results$change$significant <- outPvalue$significance
    
    # additional confidence intervals
    newresults <- unlist(resresearch(resstore, target='change$Chisq')) # should return a list
    civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$change$Chisq.ci.lower <- civals[which(names(civals) == 'lower')]
    results$change$Chisq.ci.upper <- civals[which(names(civals) == 'upper')]
    
    newresults <- unlist(resresearch(resstore, target='change$p.value')) # should return a list
    civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$change$p.value.ci.lower <- civals[which(names(civals) == 'lower')]
    if (results$change$p.value.ci.lower < 0) {results$change$p.value.ci.lower <- 0}
    results$change$p.value.ci.upper <- civals[which(names(civals) == 'upper')]
    if (results$change$p.value.ci.upper < 0) {results$change$p.value.ci.upper <- 0}
    
    newresults <- unlist(resresearch(resstore, target='change$fixed.r.squared')) # should return a list
    civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$change$fixed.r.squared.ci.lower <- civals[which(names(civals) == 'lower')]
    results$change$fixed.r.squared.ci.upper <- civals[which(names(civals) == 'upper')]
    
    newresults <- unlist(resresearch(resstore, target='change$random.r.squared')) # should return a list
    civals <- suppressWarnings(Rmisc::CI(newresults, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    results$change$random.r.squared.ci.lower <- civals[which(names(civals) == 'lower')]
    results$change$random.r.squared.ci.upper <- civals[which(names(civals) == 'upper')]
    
  }
  
  # subjects
  newresults <- resresearch(resstore, target='numparticipants') # should return a table
  results$meanofparticipants <- mean(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  results$medianofparticipants <- median(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  
  return(results)
}



resresearch <- function(resstore, target) {
  # internal function to obtain target data from list
  # returns a merged data frame if the target is a data frame
  # returns a list if the target is a list
  
  outlist <- list()
  outtable <- NULL
  booldf <- FALSE
  availablenames <- names(resstore)
  for (cAN in 1:length(availablenames)) {
    textcall <- sprintf('tempelement <- resstore$%s$%s', availablenames[cAN], target)
    eval(parse(text=textcall))
    
    if (is.data.frame(tempelement)) {
      booldf <- TRUE
      if (is.null(outtable)) {
        outtable <- tempelement
      } else {
        outtable <- rbind(outtable, tempelement)
      }
    } else {
      # store it
      textcall <- sprintf('outlist$repetition%d <- tempelement', cAN)
      eval(parse(text=textcall))
    }
  }
  if (booldf) {
    return(outtable)
  } else {
    return(outlist)
  }
}


