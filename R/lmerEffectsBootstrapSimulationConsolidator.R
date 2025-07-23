#' lmerEffectsBootstrapSimulationConsolidator
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerEffects function. This portion handles the consolidation of the models and summarizes them into a singular model.
#'
#' @param results List containing output from lmerEffects
#' @param average text parameter to indicate how the results should be collapsed. Default is using the median. Other options is mean.
#' @param reporteddata text parameter to indicate if the posthoc text reports should use actual data (actual) or should report the simulated or resampled data means.
#'
#' @return
#' \item{results}{List containing boostrapped output.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 23, 2025
#'
#' @importFrom Rmisc CI
#' @importFrom miscTools colMedians
#' @importFrom stats model.frame
#' 
#' @examples
#' \dontrun{
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
#'                covariates=NULL, planned=c('Group'), posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
#'     results <- Rmimic::lmerEffectsBootstrapSimulationCreation(results, repetitions=999, tmpdir='/tempdirectory/')
#'     results <- Rmimic::lmerEffectsBootstrapSimulationConsolidator(results, average='median', reporteddata='actual')
#'     }
#'
#' @export

lmerEffectsBootstrapSimulationConsolidator <- function(results, average='median', reporteddata='simulated') {
  
  # load data from files
  file_list <- list.files(results$futuretag$tmpdir, pattern = "^result_.*\\.RData$")
  resstore <- vector("list", length(file_list))
  for (i in seq_along(file_list)) {
    load(file.path(results$futuretag$tmpdir, file_list[i]))  # loads `result`
    resstore[[i]] <- result
  }

  Sys.sleep(1) # to make sure files are read in
  
  # summarize
  results <- resmergebootforlmerEffectsBootstrapSimulationConsolidator(results, resstore, average=average, reporteddata=reporteddata)
  
  # obtain text outputs
  if ((reporteddata == 'actual') | (reporteddata == 'raw')) {
    subtag <- 'raw'
  } else {
    if (results$futuretag$method == "resample") {
      subtag <- 'resampled'
    } else {
      subtag <- 'simulated'
    }
  }
  results <- lmerEffects2text(results, subtag=subtag, testconfidence=TRUE, significanceconfidence=TRUE)
  
  # see if posthoc adjustments are needed
  if (results$futuretag$boolposthoc) {
    if ((tolower(results$futuretag$methodofposthoccorrection) != 'false') | (tolower(results$futuretag$methodofposthoccorrection) != 'none')) {
      results <- lmerPosthocCorrection(results, method=results$futuretag$methodofposthoccorrection, studywiseAlpha=results$studywiseAlpha, FDRC=0.05)
    }
  }
  
  # tag messageout 
  if ('messageout' %in% names(results)) {
    
    outstring <- sprintf('Unstandardized effects were computed based upon bootstrapped analyses with')
    if (results$futuretag$method == "resample") {
      outstring <- sprintf('%s %d resamples', outstring, results$futuretag$repetitions)
      if (!((results$futuretag$resample_min == nrow(stats::model.frame(results$fit))) & (results$futuretag$resample_max == nrow(stats::model.frame(results$fit))))) {
        outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%)',
                             outstring, results$futuretag$resample_min, round((results$futuretag$resample_min/results$futuretag$totalsample)*100, digits=1),
                             resample_max, round((results$futuretag$resample_max/results$futuretag$totalsample)*100, digits=1))
      }
      outstring <- sprintf('%s of the original dataset (with replacement).', outstring)
    }
    if ((results$futuretag$method == "parametric") | (results$futuretag$method == "nonparametric")) {
      if (results$futuretag$method == "parametric") {
        outstring <- sprintf('%s %d datasets simulated from a multivariate normal distribution using the MASS mvrnorm function (Venables &#38; Ripley, 2002).', outstring, results$futuretag$repetitions)
      } else {
        outstring <- sprintf('%s %d datasets simulated from a multivariate non-normal distribution using the mnonr unonr function (Qu &#38; Zhang, 2020).', outstring, results$futuretag$repetitions)
      }
      outstring <- sprintf('%s For each simulation the covariance matrix was informed by', outstring)
      if (results$futuretag$subsample < 1.0) {
        outstring <- sprintf('%s a subsample of %.1f%% of the original data', outstring, round((results$futuretag$subsample)*100, digits=1))
      } else {
        outstring <- sprintf('%s the full sample of original data', outstring)
      }
      if (results$futuretag$inflation > 1.0) {
        outstring <- sprintf('%s and extrapolated to a final sample of %.0f participants', outstring, round(results$meanofparticipants, digits=0))
      } else {
        outstring <- sprintf('%s.', outstring)
      }
    }
    if (results$futuretag$method == "default") {
      outstring <- sprintf('%s %d datasets simulated by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters', outstring, repetitions)
      
      if (!((results$futuretag$resample_min == nrow(stats::model.frame(results$fit))) & (results$futuretag$resample_max == nrow(stats::model.frame(results$fit))))) {
        outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%) of the original dataset (with replacement)',
                             outstring, results$futuretag$resample_min, round((results$futuretag$resample_min/results$futuretag$totalsample)*100, digits=1),
                             results$futuretag$resample_max, round((results$futuretag$resample_max/results$futuretag$totalsample)*100, digits=1))
      }
      outstring <- sprintf('%s.', outstring)
      
      if (!is.null(results$futuretag$subsample)) {
        if (results$futuretag$subsample < 1.0) {
          outstring <- sprintf('%s For each simulation the the conditional distribution of the outcome variable was informed by', outstring)
          outstring <- sprintf('%s a subsample of %.1f%% of the original data within each between subjects factor.', outstring, round((results$futuretag$subsample)*100, digits=1))
        }
      }
      
    }
    
    results$messageout <- sprintf('%s %s', results$messageout, outstring)
  }
  
  # remove unnecessary information
  results$futuretag <- NULL
  
  return(results)
}


resmergebootforlmerEffectsBootstrapSimulationConsolidator <- function(results, resstore, average='median', reporteddata='actual') {
  
  newstats <- results$stats
  newstats$textoutput <- NA
  newstats$significance <- NA
  newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target='stats') # should return a table
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
  newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target='randomstats') # should return a table
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
  newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target='rsquared') # should return a table
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
    newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target='descriptives') # should return a table
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
  newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target='numparticipants') # should return a table
  results$meanofparticipants <- mean(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  results$medianofparticipants <- median(unlist(newresults, recursive=TRUE), na.rm=TRUE)
  
  # posthoc tests
  if ('posthoc' %in% names(results)) {
    
    uniquenames <- unique(names(results$posthoc))
    for (cUN in 1:length(uniquenames)) {
      
      # get new element
      newresults <- resresearchforlmerEffectsBootstrapSimulationConsolidator(resstore, target=sprintf('posthoc$%s', uniquenames[cUN]))
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
        newstats <- resmergebootforlmerEffectsBootstrapSimulationConsolidator(newstats, newresults, average=average, reporteddata=reporteddata)
      }
      
      # put element back
      textcall <- sprintf("results$posthoc$%s <- newstats", uniquenames[cUN])
      eval(parse(text=textcall))
      
    } # CUN
  } #has posthocs
  
  return(results)
}


resresearchforlmerEffectsBootstrapSimulationConsolidator <- function(resstore, target) {
  # function to obtain target data from list
  # returns a merged data frame if the target is a data frame
  # returns a list if the target is a list
  
  outlist <- list()
  outtable <- NULL
  booldf <- FALSE
  for (cAN in 1:length(resstore)) {
    textcall <- sprintf('tempelement <- resstore[[%d]]$%s', cAN, target)
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





