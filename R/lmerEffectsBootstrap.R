#' lmerEffectsBootstrap
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
#' @param average text parameter to indicate how the results should be collapsed. Default is using the median. Other options is mean.
#' @param reporteddata text parameter to indicate if the posthoc text reports should use actual data (actual) or should report the simulated or resampled data means.
#' @param progressbar Boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.

#' @return
#' \item{results}{List containing boostrapped output.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom dplyr slice_sample
#' @importFrom progress progress_bar
#' @importFrom Rmisc CI
#' @importFrom miscTools colMedians
#' @importFrom stats runif model.frame
#' 
#' @examples
#' \dontrun{
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
#'                covariates=NULL, planned=c('Group'), posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
#'     results <- Rmimic::lmerEffectsBootstrap(results, repetitions=999)
#'     }
#'
#' @export

lmerEffectsBootstrap <- function(results, repetitions, resample_min=NULL, resample_max=NULL, subsample=0.96, inflation=1.0, method='default', average='median', reporteddata='simulated', progressbar=TRUE) {
  
  
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
  
  
  
  
  
  # debug
  debug <- FALSE
  if (debug == TRUE) {
    data <- Rmimic::alertness
    data <- data[which(data$Condition == 'Condition1'),]
    fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = data)
    
    results <- lmerEffects(fit, subjectid = "PartID", df = "Kenward-Roger", confidenceinterval=0.95)
    results <- lmerPosthoc(results, between=c('Group'), within=c('Time', 'Condition'), covariates=c('Sex'))

    repetitions <- 10
    resample_min <- NULL
    resample_max <- NULL
    
  }
  
  totalsample <- nrow(stats::model.frame(results$fit))
  if (is.null(resample_min)) {
    resample_min <- nrow(stats::model.frame(results$fit))
  }
  if (is.null(resample_max)) {
    resample_max <- nrow(stats::model.frame(results$fit))
  }  
  boolposthoc <- FALSE
  methodofposthoccorrection <- 'false'
  if ('posthoccorrection' %in% names(results)) {
    methodofposthoccorrection <- results$posthoccorrection
  }
  if ('posthoc' %in% names(results)) {
    boolposthoc <- TRUE
    # rerun but force all posthocs
    results <- suppressWarnings(suppressMessages(lmerPosthoc(results, between=results$between, within=results$within, covariates=results$covariates, planned=results$stats$Effect, posthoclimit=results$posthoclimit, posthoccorrection='none', progressbar=FALSE)))
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
  
  # CHECK TO MAKE SURE THAT lmerEffects HAS BEEN RUN
  listofoutputsfromlmerEffects <- c('fit', 'stats', 'randomstats', 'rsquared')
  if (all(listofoutputsfromlmerEffects %in% names(results))) {
    
    if (progressbar) {
      # establish progress
      #cat(sprintf('lmerEffectsBootstrap() beginning processing '))
      countticks <- 0
      stepcounts <- 10
      checkins <- floor(seq(0, repetitions, length.out = stepcounts+1))
      pb <- progress::progress_bar$new(total = stepcounts,
                                       format = " lmerEffectsBootstrap() processing [:bar] :percent eta: :eta",
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
        newresults <- suppressWarnings(suppressMessages(lmerEffects(newfit, dependentvariable=results$dependentvariable, subjectid=results$subjectid, within=results$within, df = results$df, confidenceinterval=results$confidenceinterval, studywiseAlpha=results$studywiseAlpha, suppresstext=TRUE, smp=smp)))
        # compute posthoc if previously run - function should have stored the necessary information
        if (boolposthoc) {
          newresults <- invisible(suppressWarnings(suppressMessages(lmerPosthoc(newresults, between=results$between, within=results$within, covariates=results$covariates, planned=results$stats$Effect, posthoclimit=results$posthoclimit, calltype='subprocess', posthoccorrection='none'))))
        }
        newresults$descriptives <- lmerEffects_simpledesc(newresults) # store these
        
        # remove model
        newresults$fit <- NULL
        
        # store it
        textcall <- sprintf('resstore$repetition%d <- newresults', cN)
        eval(parse(text=textcall))
        
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
    results <- tryCatch({
      results <- resmergeboot(results, resstore, average=average, reporteddata=reporteddata)
    }, error = function(e) {
      cat(sprintf('lmerEffectsBootstrap - summarize failure\n'))
      results <- NULL
    })
    
    # obtain text outputs
    if ((reporteddata == 'actual') | (reporteddata == 'raw')) {
      subtag <- 'raw'
    } else {
      if (method == "resample") {
        subtag <- 'resampled'
      } else {
        subtag <- 'simulated'
      }
    }
    
    results <- tryCatch({
      results <- lmerEffects2text(results, subtag=subtag, testconfidence=TRUE, significanceconfidence=TRUE)
    }, error = function(e) {
      cat(sprintf('lmerEffectsBootstrap - lmerEffects2text failure\n'))
      results <- NULL
    })
    
    # see if posthoc adjustments are needed
    if (boolposthoc) {
      if ((tolower(methodofposthoccorrection) != 'false') | (tolower(methodofposthoccorrection) != 'none')) {
        results <- lmerPosthocCorrection(results, method=methodofposthoccorrection, studywiseAlpha=results$studywiseAlpha, FDRC=0.05)
      }
    }
    
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


resmergeboot <- function(results, resstore, average='median', reporteddata='actual') {
  
  newstats <- results$stats
  newstats$textoutput <- NA
  newstats$significance <- NA
  newresults <- resresearch(resstore, target='stats') # should return a table
  # fixed effects
  for (cR in 1:nrow(newstats)) {
    tempdbs <- newresults[which(newresults$Effect == newstats$Effect[cR]),]
    colsofinterest <- c("DFn", "DFd", "SSn", "SSd", "SSe", "F.value", "p.value", "partialetasquared", "fsquared","fsquared.ci.lower", "fsquared.ci.upper")
    for (cC in 1:length(colsofinterest)) {
      if (average == 'median') {
        newstats[cR,colsofinterest[cC]] <- median(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      } else {
        newstats[cR,colsofinterest[cC]] <- mean(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      }
    }
    
    # significance test
    outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
    newstats$significance[cR] <- outPvalue$significance
    
    # additional confidence intervals
    civals <- suppressWarnings(Rmisc::CI(tempdbs$F.value, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    newstats$F.value.ci.lower[cR] <- civals[which(names(civals) == 'lower')]
    newstats$F.value.ci.upper[cR] <- civals[which(names(civals) == 'upper')]
    
    civals <- suppressWarnings(Rmisc::CI(tempdbs$p.value, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    newstats$p.value.ci.lower[cR] <- civals[which(names(civals) == 'lower')]
    if (newstats$p.value.ci.lower[cR] < 0) {newstats$p.value.ci.lower[cR] <- 0}
    newstats$p.value.ci.upper[cR] <- civals[which(names(civals) == 'upper')]
    if (newstats$p.value.ci.upper[cR] < 0) {newstats$p.value.ci.upper[cR] <- 0}
  }
  results$stats <- newstats # swap out
  
  # random effects
  newstats <- results$randomstats
  newstats$textoutput <- NA
  newstats$significance <- NA
  newresults <- resresearch(resstore, target='randomstats') # should return a table
  # create empty variables
  newstats$p.value.ci.lower <- NA
  newstats$p.value.ci.upper <- NA
  for (cR in 1:nrow(newstats)) {
    tempdbs <- newresults[which(newresults$Effect == newstats$Effect[cR]),]
    colsofinterest <- c("DF", "LogLikelihood", "LRT", "p.value")
    for (cC in 1:length(colsofinterest)) {
      if (average == 'median') {
        newstats[cR,colsofinterest[cC]] <- median(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      } else {
        newstats[cR,colsofinterest[cC]] <- mean(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      }
    }
    # significance test
    outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
    newstats$significance[cR] <- outPvalue$significance
    
    civals <- suppressWarnings(Rmisc::CI(tempdbs$p.value, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
    newstats$p.value.ci.lower[cR] <- civals[which(names(civals) == 'lower')]
    if (newstats$p.value.ci.lower[cR] < 0) {newstats$p.value.ci.lower[cR] <- 0}
    newstats$p.value.ci.upper[cR] <- civals[which(names(civals) == 'upper')]
    if (newstats$p.value.ci.upper[cR] < 0) {newstats$p.value.ci.upper[cR] <- 0}
  }
  results$randomstats <- newstats
  
  # rsquared
  newstats <- results$rsquared
  newresults <- resresearch(resstore, target='rsquared') # should return a table
  for (cR in 1:nrow(newstats)) {
    tempdbs <- newresults[which(newresults$portion == newstats$portion[cR]),]
    colsofinterest <- c("effects", "ci.lower", "ci.upper")
    for (cC in 1:length(colsofinterest)) {
      if (average == 'median') {
        newstats[cR,colsofinterest[cC]] <- median(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      } else {
        newstats[cR,colsofinterest[cC]] <- mean(tempdbs[,colsofinterest[cC]], na.rm=TRUE)
      }
    }
  }
  results$rsquared <- newstats
  
  # descriptives
  if ('descriptives' %in% names(results)) {
    newstats <- results$descriptives
    newstats$DistributionData <- NA
    newresults <- resresearch(resstore, target='descriptives') # should return a table
    for (cR in 1:nrow(newstats)) {
      tempdbs <- newresults[which(newresults$Group == newstats$Group[cR]),]
      colsofinterest <- c("N", "Missing", "Mean", "Median", "SD", "SE")
      for (cC in 1:length(colsofinterest)) {
        if (average == 'median') {
          newstats[cR,colsofinterest[cC]] <- median(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
        } else {
          newstats[cR,colsofinterest[cC]] <- mean(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
        }
      }
      distributiontable <- data.frame(matrix(NA, nrow=nrow(tempdbs), ncol=64))
      for (csR in 1:nrow(tempdbs)) {
        textcall <- sprintf('distributiontable[csR,] <- c(%s)', tempdbs$DistributionData[csR])
        eval(parse(text=textcall))
      }
      for (cC in 1:ncol(distributiontable)) {
        distributiontable[,cC] <- as.numeric(distributiontable[,cC])
      }
      if (average == 'median') {
        newstats$DistributionData[cR] <- paste(miscTools::colMedians(distributiontable, na.rm=TRUE), collapse = ",")
      } else {
        newstats$DistributionData[cR] <- paste(colMeans(distributiontable, na.rm=TRUE), collapse = ",")
      }
      decisionindx <- which(tempdbs$DistributionDecision == "Normal")
      if (length(decisionindx) > 0) {
        if ((length(decisionindx) / nrow(tempdbs)) >= 0.6) {
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
  }
  
  # subjects
  newresults <- resresearch(resstore, target='numparticipants') # should return a table
  results$meanofparticipants <- mean(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  results$medianofparticipants <- median(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  
  # posthoc tests
  if ('posthoc' %in% names(results)) {
    
    uniquenames <- unique(names(results$posthoc))
    for (cUN in 1:length(uniquenames)) {
      
      # get new element
      newresults <- resresearch(resstore, target=sprintf('posthoc$%s', uniquenames[cUN]))
      textcall <- sprintf("newstats <- results$posthoc$%s", uniquenames[cUN])
      eval(parse(text=textcall))
      
      if (is.data.frame(newresults)) {
        # element is a table
        newstats$textoutput <- NA
        newstats$significant <- NA
        for (cR in 1:nrow(newstats)) {
          tempdbs <- newresults[which(newresults$idtag == newstats$idtag[cR] & newresults$hold == newstats$hold[cR]),]
          if (nrow(tempdbs) == 0) {
            # data simulation likely swapped position of c1 and c2
            tempvect <- stringr::str_split(newstats$idtag[cR], "-")[[1]]
            newtempvect <- tempvect
            newtempvect[length(tempvect)] <- tempvect[length(tempvect)-2]
            newtempvect[length(tempvect)-2] <- tempvect[length(tempvect)]
            newstats$idtag[cR] <- paste0(newtempvect, collapse="-")
            c1name <- newstats$C1name
            c2name <- newstats$C2name
            newstats$C1name <- c2name
            newstats$C2name <- c1name
            tempdbs <- newresults[which(newresults$idtag == newstats$idtag[cR] & newresults$hold == newstats$hold[cR]),]
          }
          colsofinterest <- c("df", "t.ratio", "p.value", "effectsize", "effectsize.conf.int.lower", "effectsize.conf.int.upper", "correlation")
          for (cC in 1:length(colsofinterest)) {
            if (average == 'median') {
              newstats[cR,colsofinterest[cC]] <- median(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
            } else {
              newstats[cR,colsofinterest[cC]] <- mean(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
            }
          }
          
          # significance test
          outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
          newstats$significant[cR] <- outPvalue$significance
          
          # additional confidence intervals
          tempvectone <- as.numeric(tempdbs$t.ratio)
          tempvectone <- tempvectone[which(!is.na(tempvectone))]
          if (length(tempvectone) > 0) {
            civals <- suppressWarnings(Rmisc::CI(tempvectone, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
            newstats$t.conf.int.lower[cR] <- civals[which(names(civals) == 'lower')]
            newstats$t.conf.int.upper[cR] <- civals[which(names(civals) == 'upper')]
            if (!is.na(newstats$t.conf.int.lower[cR])) {
              if (newstats$t.conf.int.lower[cR] > newstats$t.ratio[cR]) {
                newstats$t.conf.int.lower[cR] <- newstats$t.ratio[cR]
              }
            }
            if (!is.na(newstats$t.conf.int.upper[cR])) {
              if (newstats$t.conf.int.upper[cR] < newstats$t.ratio[cR]) {
                newstats$t.conf.int.upper[cR] <- newstats$t.ratio[cR]
              }
            }
          } else {
            newstats$t.conf.int.lower[cR] <- NA
            newstats$t.conf.int.upper[cR] <- NA
          }
          
          tempvectone <- as.numeric(tempdbs$p.value)
          tempvectone <- tempvectone[which(!is.na(tempvectone))]
          if (length(tempvectone) > 0) {
            civals <- suppressWarnings(Rmisc::CI(tempvectone, ci = (results$confidenceinterval-results$studywiseAlpha))) # one sided
            newstats$p.conf.int.lower[cR] <- civals[which(names(civals) == 'lower')]
            if (!is.na(newstats$p.conf.int.lower[cR])) {
              if (newstats$p.conf.int.lower[cR] < 0) {newstats$p.conf.int.lower[cR] <- 0}
            }
            newstats$p.conf.int.upper[cR] <- civals[which(names(civals) == 'upper')]
            if (!is.na(newstats$p.conf.int.upper[cR])) {
              if (newstats$p.conf.int.upper[cR] < 0) {newstats$p.conf.int.upper[cR] <- 0}
            }
            if (!is.na(newstats$p.conf.int.lower[cR])) {
              if (newstats$p.conf.int.lower[cR] > newstats$p.value[cR]) {
                newstats$p.conf.int.lower[cR] <- newstats$p.value[cR]
              }
            }
            if (!is.na(newstats$p.conf.int.upper[cR])) {
              if (newstats$p.conf.int.upper[cR] < newstats$p.value[cR]) {
                newstats$p.conf.int.upper[cR] <- newstats$p.value[cR]
              }
            }
          } else {
            newstats$p.conf.int.lower[cR] <- NA
            newstats$p.conf.int.upper[cR] <- NA
          }
          
          if (reporteddata != 'actual') {
            colsofinterest <- c("C1n", "C1mean", "C1sd", "C2n", "C2mean", "C2sd")
            for (cC in 1:length(colsofinterest)) {
              if (average == 'median') {
                newstats[cR,colsofinterest[cC]] <- median(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
              } else {
                newstats[cR,colsofinterest[cC]] <- mean(as.numeric(tempdbs[,colsofinterest[cC]]), na.rm=TRUE)
              }
            }  
          }
        }
        
      } else {
        
        # recursive call
        newstats <- resmergeboot(newstats, newresults, average=average, reporteddata=reporteddata)
        
      }
      
      # put element back
      textcall <- sprintf("results$posthoc$%s <- newstats", uniquenames[cUN])
      eval(parse(text=textcall))
      
    } # CUN
  } #has posthocs
  
  return(results)
}






