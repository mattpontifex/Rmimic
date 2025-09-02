#' lmerPosthoc
#'
#' @description Posthoc decomposition of univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. Interactions are decomposed multiple ways (A holding B, B holding A) and superseeding interactions suppress lower level effects tests (no posthoc test of A:B if A:B:C is significant). Tests of A:B can still be obtained using the planned parameter if desired.
#'
#' @param results list output from the lmerEffects function
#' @param between text list specifying the between subjects variables
#' @param within text list specifying any variables that should be considered as within subjects. If null assumes all variables are between subjects.
#' @param covariates text list specifying the variables that should be considered as non interactive covariates. 
#' @param planned text list specifying any effect to show the post-hoc comparisons even if they are not significant.
#' @param suppress text list specifying any effect to ignore even if it is significant.
#' @param posthoccorrection text parameter to indicate what post-hoc comparisons should be performed. Default is False Discovery Rate Control. Other options are Bonferroni, Holm-Bonferroni, Sidak. None will skip post hoc corrections.
#' @param bootstrap list paramater to indicate if bootstrapping should be performed. Default is NULL. List should provide the call elements for lmerEffectsBootstrapANOVA and be in the form list('repetitions' = 10000, 'subsample'=0.96, 'inflation'=1.0, 'method'='default')
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param posthoclimit integer specifying how many factors in an interaction should be broken down. Default is 6 indicating posthoc results will be provided for up to a 6 way interaction.
#' @param progressbar boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.
#'
#' @return results list output from the lmerEffects function
#' \item{posthoc}{list output containing the posthoc test results}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#' 
#' @importFrom stats model.frame
#' 
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerPosthoc(results, dependentvariable = "Alertness", subjectid = "PartID", 
#'                between=c('Group'), within=c('Time'), covariates=NULL,
#'                planned=c('Group'), df = "Kenward-Roger",
#'                posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
#'                
#' @export

lmerPosthoc <- function(results, between=NULL, within=NULL, covariates=NULL, dependentvariable=NULL, subjectid=NULL, df = NULL, planned=NULL, suppress=NULL, posthoccorrection=NULL, bootstrap=NULL, confidenceinterval=NULL, studywiseAlpha=NULL, posthoclimit=6, calltype=NULL, verbose=FALSE, progressbar=TRUE, ...) {
  
  debug <- FALSE
  
  if (debug) {
    within=NULL
    between=NULL
    covariates=NULL
    confidenceinterval=NULL
    studywiseAlpha=NULL
    
    within=c('Time')
    between=c('Group')
    covariates=c('Sex')
    
    data <- Rmimic::alertness
    data <- data[which(data$Condition == 'Condition2'),]
    #data <- data[which(data$Time != 'Time2'),]
    fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = data)
    
    results <- lmerEffects(fit, subjectid = "PartID", df = "Kenward-Roger", confidenceinterval=0.95)
    
    #fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = data)
    
    #data <- Rmimic::alertness
    #data <- data[which(data$Condition == 'Condition2'),]
    #data$NewGroup <- NA
    #data$NewGroup[which(data$Group == 'Group1')] <- 1
    #data$NewGroup[which(data$Group == 'Group2')] <- 2
    #data$Group <- as.numeric(data$NewGroup)
    #fit <- lmerTest::lmer(Alertness ~ Group*Time + Sex + (1 | PartID), data = data)
    
    #results <- lmerEffects(fit, subjectid = "PartID", df = "Kenward-Roger", confidenceinterval=0.95)
    
  }
  
  # CHECK TO MAKE SURE THAT lmerEffects HAS BEEN RUN
  listofoutputsfromlmerEffects <- c('fit', 'stats', 'randomstats', 'rsquared')
  if (!all(listofoutputsfromlmerEffects %in% names(results))) {
    if (!is.null(bootstrap)) {
      df = "Shattertwaite" # computation time becomes extreme otherwise
    }
    starttime <- Sys.time()
    # try running lmerEffects
    results <- tryCatch({
      results <- invisible(suppressWarnings(suppressMessages(lmerEffects(results, dependentvariable=dependentvariable, subjectid=subjectid, within=within, df=df, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, suppresstext=FALSE, smp=NULL, verbose=verbose))))
    }, error = function(e) {
      results <- NULL
    })
    if (verbose) {
      cat(sprintf('  lmerPosthoc(): time to run lmerEffects - %.2f sec\n', Sys.time() - starttime))
    }
  }
  
  if (is.null(studywiseAlpha)) {
    studywiseAlpha <- tryCatch({
      studywiseAlpha <- results$studywiseAlpha
    }, error = function(e) {
      studywiseAlpha <- NULL
    })
    if (is.null(studywiseAlpha)) {
      studywiseAlpha <- 0.05
      results$studywiseAlpha <- studywiseAlpha # store it
    }
  }
  if (is.null(confidenceinterval)) {
    confidenceinterval <- tryCatch({
      confidenceinterval <- results$confidenceinterval
    }, error = function(e) {
      confidenceinterval <- NULL
    })
    if (is.null(confidenceinterval)) {
      confidenceinterval <- 0.95
      results$confidenceinterval <- confidenceinterval # store it
    }
  }
  
  if (!is.null(bootstrap)) {
    if (any(bootstrap)) {
      df = "Shattertwaite" # computation time becomes extreme otherwise
    }
  }
  if (is.null(df)) {
    df <- tryCatch({
      df <- results$df
    }, error = function(e) {
      df <- NULL
    })
    if (is.null(df)) {
      df <- "Kenward-Roger"
      results$df <- df # store it
    }
  }
  
  if (is.null(posthoccorrection)) {
    posthoccorrection <- tryCatch({
      posthoccorrection <- results$posthoccorrection
    }, error = function(e) {
      posthoccorrection <- NULL
    })
    if (is.null(posthoccorrection)) {
      posthoccorrection <- "False Discovery Rate Control"
    }
  }
  results$posthoccorrection <- posthoccorrection # store it
  
  if (!is.null(df)) {
    if (toupper(df) == toupper("Kenward-Roger")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("KR")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("S")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("satterthwaite")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("Traditional")) {
      df = "Traditional"
    }
  } else {
    df = "Kenward-Roger"
  }
  emmeansdf <- tolower(df)
  if (df == "Kenward-Roger") {
    emmeansdf <- 'kenward-roger'
  } else {
    emmeansdf <- 'satterthwaite'
  }
  
  tempdbs <- tryCatch({
    tempdbs <- stats::model.frame(results$fit)
  }, error = function(e) {
    tempdbs <- NULL
  })
  if (!is.null(tempdbs)) {
    
    if (is.null(dependentvariable)) {
      dependentvariable <- tryCatch({
        dependentvariable <- results$dependentvariable
      }, error = function(e) {
        dependentvariable <- NULL
      })
      if (is.null(dependentvariable)) {
        dependentvariable <- colnames(tempdbs)[1] # assume standard format
        results$dependentvariable <- dependentvariable # store it
      }
    }
    
    if (is.null(subjectid)) {
      subjectid <- tryCatch({
        subjectid <- results$subjectid
      }, error = function(e) {
        if (!is.null(within)) {
          warning('lmerPosthoc(): Within subjects factor specified but no factor specified for Subject.')
        }
        subjectid <- NULL
      })
    }
    
    # store some information
    withintest <- tryCatch({
      withintest <- results$within
    }, error = function(e) {
      withintest <- NULL
    })
    if (is.null(withintest)) {
      results$within <- within
    }
    covariatestest <- tryCatch({
      covariatestest <- results$covariates
    }, error = function(e) {
      covariatestest <- NULL
    })
    if (is.null(covariatestest)) {
      results$covariates <- covariates
    }
    betweentest <- tryCatch({
      betweentest <- results$between
    }, error = function(e) {
      betweentest <- NULL
    })
    if (is.null(betweentest)) {
      results$between <- between
    }
    posthoclimittest <- tryCatch({
      posthoclimittest <- results$posthoclimit
    }, error = function(e) {
      posthoclimittest <- NULL
    })
    if (is.null(posthoclimittest)) {
      results$posthoclimit <- posthoclimit
    }
    plannedtest <- tryCatch({
      plannedtest <- results$planned
    }, error = function(e) {
      plannedtest <- NULL
    })
    if (is.null(plannedtest)) {
      results$planned <- planned
    }
    suppresstest <- tryCatch({
      suppresstest <- results$suppress
    }, error = function(e) {
      suppresstest <- NULL
    })
    if (is.null(suppresstest)) {
      results$suppress <- suppress
    }
    
    # check to see if bootstrapping should be performed
    if (!is.null(bootstrap)) {
      results <- tryCatch({
        if (any(bootstrap)) {
          results <- lmerEffectsBootstrapANOVA(results, bootstrap)
        } else {
          results <- results
        }
      }, error = function(e) {
        results <- results
      })
    }
    
    workingdbs <- tryCatch({
      workingdbs <- results$stats
    }, error = function(e) {
      workingdbs <- NULL
    })
    if (!is.null(workingdbs)) {
      
      if (!('significance' %in% colnames(workingdbs))) {
        workingdbs$significance <- NA
        for (cRsig in 1:nrow(workingdbs)) {
          workingdbs$significance[cRsig] <- fuzzyP(workingdbs$p.value[cRsig], studywiseAlpha)$significance
        }
      }
      
      workingdbs$decompose <- TRUE # set all to decompose
      workingdbs$decompose <- workingdbs$significance # only significant show
      superseedinteraction <- tryCatch({
        
        # look for superseeding interactions
        if (length(which(workingdbs$decompose)) > 0) {
          uniquedecomps <- workingdbs$Effect[which(workingdbs$decompose)]
          for (cUDC in 1:length(uniquedecomps)) {
            for (cUDC2 in 1:length(uniquedecomps)) {
              comp1 <- uniquedecomps[cUDC]
              comp2 <- uniquedecomps[cUDC2]
              if (comp1 != comp2) {
                factorsinvolved1 <- unlist(strsplit(as.character(comp1),"[:]"))
                factorsinvolved2 <- unlist(strsplit(as.character(comp2),"[:]"))
                if (all(factorsinvolved1 %in% factorsinvolved2)) {
                  # all of factor 1 is covered by factor 2
                  workingdbs$decompose[which(workingdbs$Effect == comp1)] <- FALSE # turn off decomposition
                }
              }
            }
          }
        }
        
        superseedinteraction <- 1
      }, error = function(e) {
        superseedinteraction <- NULL
      })
      
      # do some quick data checking
      for (currentAnovaLine in 1:nrow(workingdbs)) {
        
        # see if the effect is with a covariate
        factorsinvolved <- unlist(strsplit(as.character(workingdbs$Effect[currentAnovaLine]),"[:]"))
        if (!is.null(covariates)) {
          if (any(covariates %in% factorsinvolved)) {
            # if the interaction is with a covariate then skip
            if (length(factorsinvolved) > 1) {
              workingdbs$decompose[currentAnovaLine] <- FALSE
            }
          }
        }
        
        # see if it was requested to be suppressed
        if (workingdbs$Effect[currentAnovaLine] %in% suppress) {
          workingdbs$decompose[currentAnovaLine] <- FALSE
        }
        
        # see if it is planned
        if (workingdbs$Effect[currentAnovaLine] %in% planned) {
          workingdbs$decompose[currentAnovaLine] <- TRUE
        }
        
        # see if it is a simple decomp
        if (length(factorsinvolved) == 1) {
          tempdbs <- stats::model.frame(results$fit)
          if (length(unique(tempdbs[,factorsinvolved])) < 3) {
            workingdbs$decompose[currentAnovaLine] <- TRUE
          }
        }
      }
      
      if (!is.null(calltype)) {
        if (!any(workingdbs$decompose)) {
          # subprocess but nothing is significant
          workingdbs$decompose[which(workingdbs$factorsinvolved == 1)] <- TRUE
        }
      }
      
      for (currentAnovaLine in 1:nrow(workingdbs)) {
        
        subplanned <- NULL
        if (workingdbs$Effect[currentAnovaLine] %in% planned) {
          subplanned <- workingdbs$Effect # if planned all decompositions go in for the subprocess
        }
        
        if (workingdbs$decompose[currentAnovaLine]) {
          starttime <- Sys.time()
          
          # obtain breakdowns
          tempresult <- results
          tempresult$planned <- subplanned
          tempresult$covariates <- covariates
          tempresult$progressbar <- progressbar
          tempresult$bootstrap <- bootstrap
          tempresult <- invisible(suppressWarnings(suppressMessages(lmerPosthocsubprocess(tempresult, workingdbs$Effect[currentAnovaLine]))))
          
          if (verbose) {
            cat(sprintf('  lmerPosthoc(): time to decompose effect %d - %.2f sec\n', currentAnovaLine, Sys.time() - starttime))
          }
          
          if (length(names(tempresult)) > 0) {
            # store results
            if ('posthoc' %in% names(results)) {
              results$posthoc[names(tempresult)] <- tempresult
            } else {
              results$posthoc <- list()
              results$posthoc[names(tempresult)] <- tempresult
            }
          }
        }
        
      } # anova loop
      
    } # stats extraction
  } # fit exists
  
  # obtain text outputs
  starttime <- Sys.time()
  if (is.null(bootstrap)) {
    results <- lmerEffects2text(results)
  } else {
    subtag <- 'simulated'
    if ("reporteddata" %in% names(bootstrap)) {
      if ((bootstrap$reporteddata == 'actual') | (bootstrap$reporteddata == 'raw')) {
        subtag <- 'raw'
      } else if (bootstrap$reporteddata == "resample") {
        subtag <- 'resampled'
      }
    }
    results <- lmerEffects2text(results, subtag=subtag, testconfidence=TRUE, significanceconfidence=TRUE)
  }
  
  if (verbose) {
    cat(sprintf('  lmerPosthoc(): time to process text - %.2f sec\n', Sys.time() - starttime))
  }
  
  if (is.null(calltype)) {
    # perform posthoc correction
    if (!is.null(posthoccorrection)) {
      if ((tolower(posthoccorrection) != 'false') | (tolower(posthoccorrection) != 'none')) {
        
        starttime <- Sys.time()
        results <- lmerPosthocCorrection(results, method=posthoccorrection, studywiseAlpha=studywiseAlpha, FDRC=0.05)
        if (verbose) {
          cat(sprintf('  lmerPosthoc(): time to compute post hoc corrections - %.2f sec\n', Sys.time() - starttime))
        }
      }
    }
  }
  
  return(results) 
}