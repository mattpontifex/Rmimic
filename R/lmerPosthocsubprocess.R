#' lmerPosthocsubprocess
#'
#' @description Subfunction of Rmimic::lmerPosthoc. Performs posthoc tests.
#'
#' @param fit modelfit from the lmer function
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label.
#' @param effectofinterest text specifying the effect label to decompose
#' @param within text list specifying any variables that should be considered as within subjects
#' @param between text list specifying the between subjects variables
#' @param covariates text list specifying the variables that should be considered as non interactive covariates. 
#' @param planned text list specifying any effect to show the post-hoc comparisons even if they are not significant.
#' @param df Parameter to indicate what degrees of freedom approximation should be used. Default is Kenward-Roger. Other options are Shattertwaite or Traditional.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param posthoclimit integer specifying how many factors in an interaction should be broken down. Default is 6 indicating posthoc results will be provided for up to a 6 way interaction.
#' @param progressbar boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 22, 2025
#'
#' @importFrom progress progress_bar
#' @importFrom stringr str_split str_remove_all
#' @importFrom pkgcond suppress_conditions
#' @importFrom lmerTest lmer
#' @importFrom emmeans emmeans
#' @importFrom doBy summaryBy
#' @importFrom MBESS conf.limits.nct
#' @importFrom stats cor.test formula model.frame
#' @importFrom WRS2 pbcor
#' 
#' @export

lmerPosthocsubprocess <- function(fit, dependentvariable, subjectid, effectofinterest, within, between=NULL, covariates=NULL, planned=NULL, df=NULL, confidenceinterval=0.95, studywiseAlpha=0.05, posthoclimit=6, progressbar=TRUE) {

  debug <- FALSE
  if (debug) {
  
    dependentvariable <- 'Alertness'
    subjectid <- 'PartID'
    effectofinterest <- "Group:Time"
    within <- c('Time')
    df=NULL
    confidenceinterval=0.95
    studywiseAlpha=0.05
  
  }
  
  if (!is.null(df)) {
    if (toupper(df) == toupper("Kenward-Roger")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("Traditional")) {
      df = "Traditional"
    }
  } else {
    df = "Kenward-Roger"
  }
  emmeansdf <- tolower(df)
  if (emmeansdf != "satterthwaite") {
    emmeansdf <- 'kenward-roger'
  }
  
  # check that we are dealing with a lmer fit
  randomformula <- tryCatch({
    randomformula <- Reduce(paste, deparse(formula(fit, random.only = TRUE)))
  }, error = function(e) {
    randomformula <- NULL
  })
  if (!is.null(randomformula)) {
    # get data
    tempdbs <- tryCatch({
      tempdbs <- stats::model.frame(fit)
    }, error = function(e) {
      tempdbs <- NULL
    })
    
    factorsinvolved <- unlist(strsplit(as.character(effectofinterest),"[:]"))
    factorsinvolvedL <- length(factorsinvolved)
    
    # make sure this is something that we can do
    if (factorsinvolvedL < posthoclimit) {
        
      # see what we are working with
      factorlengthmatrix <- data.frame(matrix(NA,nrow=1, ncol=factorsinvolvedL))
      colnames(factorlengthmatrix) <- factorsinvolved
      for (cB in 1:factorsinvolvedL) {
        factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(tempdbs[,factorsinvolved[cB]]))))
      }
      
      res <- list()
      
      if (progressbar) {
        # establish progress
        formtext <- sprintf(" lmerPosthoc() decomposing %s [:bar] :percent eta: :eta", effectofinterest)
        #cat(sprintf('lmerPosthoc() beginning processing '))
        countticks <- 0
        stepcounts <- factorsinvolvedL
        pb <- progress::progress_bar$new(total = stepcounts,
                                         format = formtext,
                                         clear = TRUE, width= 120)
        pb$tick(0)
      }
      
      # loop through each factor
      for (currentFactorinvolved in 1:factorsinvolvedL) {
        currentfactor <- colnames(factorlengthmatrix)[currentFactorinvolved]
        currentfactorlevelsinvolved <- unique(unlist(as.character(tempdbs[,currentfactor])))
        otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactor)]
        factortag <- currentfactor
        
        if (factorsinvolvedL > 1) {
          if (length(otherfactorsinvolved) == 1) {
            decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactor[1])
            factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), otherfactorsinvolved[1], currentfactor[1])
          } else {
            decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactor[1])
            factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), 
                                 paste(otherfactorsinvolved, collapse=sprintf("_by_")), currentfactor[1])
          }
        } else {
          decomptext <- ''
        }
        decompdir <- paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 "))
        
        # establish output
        colsofinterest <- c("decomp", "hold","contrast", "df",
                            "t.ratio", 't.conf.int.lower', 't.conf.int.upper',
                            "p.value", 'p.conf.int.lower', 'p.conf.int.upper',
                            'effectsize', 'effectsize.conf.int.lower', 'effectsize.conf.int.upper', 'significant', 'correlation',
                            'C1name', 'C1n', 'C1mean', 'C1sd',
                            'C2name', 'C2n', 'C2mean', 'C2sd', 'idtag')
        outtable <- data.frame(matrix(NA, nrow=0, ncol=length(colsofinterest)))
        colnames(outtable) <- colsofinterest
        
        # check if anova is required
        requireanova <- FALSE
        if (factorsinvolvedL > 2) {
          requireanova <- TRUE # 3+ way interaction
        } else {
         # one or two way
         if (length(otherfactorsinvolved) > 0) {
           if (factorlengthmatrix[1,otherfactorsinvolved] > 2) {
             requireanova <- TRUE # more than 3 levels
           }
         }
        }
        
        if (requireanova) {  
          # ANOVA required
          for (cD in 1:length(currentfactorlevelsinvolved)) {
            # subset data
            subworkingdatabase <- tempdbs
            subworkingdatabase <- subworkingdatabase[which(subworkingdatabase[,currentfactor[1]] == currentfactorlevelsinvolved[cD]),]
            
            currentfactorlevelstring <- currentfactorlevelsinvolved[cD]
            # see if factor is just a number
            if (!grepl("\\D", currentfactorlevelstring)) {
              currentfactorlevelstring <- sprintf('%s%s', currentfactor[1], currentfactorlevelstring)
            }
            decompconst <- sprintf('When %s is %s', currentfactor[1], currentfactorlevelstring)
            
            # validate if random factors can be included
            randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
            randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]][2]
            randomformula <- stringr::str_split(randomformula, '[+]')[[1]]
            for (cRF in 1:length(randomformula)) {
              # This regex captures words composed of letters, digits, and underscores
              variables_in_formula <- unlist(regmatches(randomformula[cRF], gregexpr("\\b\\w+\\b", randomformula[cRF])))
              tempvect <- intersect(variables_in_formula, colnames(subworkingdatabase))
              boolpass <- rep_len(TRUE, length(tempvect))
              for (cRFtemp in 1:length(tempvect)) {
                if (length(unique(subworkingdatabase[,tempvect[cRFtemp]])) < 2) {
                  boolpass[cRFtemp] <- FALSE
                }
              }
              if (!all(boolpass)) {
                # insufficient levels for this random factor to remain
                boolresolved <- FALSE
                
                # can you put a constant in but keep the factor
                potentialswitch <- stringr::str_replace_all(randomformula[cRF], tempvect[which(!boolpass)], '1')
                potentialswitch <- stringr::str_replace_all(potentialswitch, ':1', '')
                potentialswitch <- stringr::str_replace_all(potentialswitch, '1:', '')
                if (!(potentialswitch %in% randomformula)) {
                  if (!(potentialswitch == '(1|1)')) {
                    randomformula[cRF] <- potentialswitch
                    boolresolved <- TRUE
                  }
                }
                # if not resolved just remove it
                if (!boolresolved) {
                  randomformula[cRF] <- ""
                }
              }
            }
            randomformula <- randomformula[which(randomformula != '')]
            randomformula <- paste(randomformula, collapse='+')
            
            if (factorsinvolvedL > 1) {
              fixedformula <- paste(otherfactorsinvolved, collapse=sprintf("*")) # since this is an interaction decomposition
              if (length(otherfactorsinvolved) == 1) {
                #decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), otherfactorsinvolved[1], currentfactorlevelstring)
              } else {
                #decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), 
                                     paste(otherfactorsinvolved, collapse=sprintf("_by_")), currentfactorlevelstring)
              }
            } else {
              fixedformula <- currentfactor[1]
              factortag <- currentfactor
            }
              
            # run the anova on the factor of interest using the subsetted data
            smp <- subworkingdatabase
            subfit <- tryCatch({
              textcall <- sprintf('subfit <- pkgcond::suppress_conditions(lmerTest::lmer(formula = %s ~ %s + %s, data = smp))', dependentvariable[1], fixedformula, randomformula)
              subfit <- eval(parse(text=textcall))
            }, error = function(e) {
              cat(sprintf('\nThis error is usually the result of too many random factors specified for the interaction.\nUnable to run decomposition anova with the formula: %s ~ %s + %s\n', dependentvariable[1], fixedformula, randomformula))
              subfit <- NULL
            })
            if (!is.null(subfit)) {
              # get table
              subresults <- lmerEffects(subfit, dependentvariable=dependentvariable, subjectid = subjectid, within=within, df = df, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, suppresstext=FALSE, smp=smp)
              # get posthoc
              subresults <- lmerPosthoc(subresults, between=between, within=within, covariates=covariates, planned=planned, posthoccorrection='none', confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, posthoclimit=posthoclimit-1, calltype='subprocess', progressbar=progressbar)
              
              # need to output subresults to res
              if (length(names(subresults)) > 0) {
                textcall <- sprintf("res$ANOVA_%s <- subresults", factortag)
                eval(parse(text=textcall))
              }
            } # submodel
          } # each factor level involved
        } else {
          # compute all pairwise contrasts
          posthoctemptestfull <- tryCatch({
            hold <- ''
            if (factorsinvolvedL > 1) {
              hold <- sprintf("by = '%s', ", currentfactor)
            }
            textcall <- sprintf("posthoctemptestfull <- pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fit, ~ %s, adjust = 'none', mode='%s'), %sadjust='none')))", paste(effectofinterest, collapse=sprintf("*")), emmeansdf, hold)
            eval(parse(text=textcall))
            posthoctemptestfull <- as.data.frame(posthoctemptestfull)
          }, error = function(e) {
            posthoctemptestfull <- NULL
          })
          if (!is.null(posthoctemptestfull)) {
            
            #"The observed noncentrality parameter of the noncentral t-distribution has exceeded 37.62
            # in magnitude (R's limitation for accurate probabilities from the noncentral t-distribution)
            # in the function's iterative search for the appropriate value(s). The results may be fine, 
            # but they might be inaccurate; use caution."
            if ('t.ratio' %in% names(posthoctemptestfull)) {
              posthoctemptestfull$t.ratio <- as.numeric(posthoctemptestfull$t.ratio)
            } else {
              # formula is the same, the p value is looked up on a different table due to asymptoc degrees of freedom
              posthoctemptestfull$t.ratio <- as.numeric(posthoctemptestfull$z.ratio)
            }
            posthoctemptestfull$t.ratio[which(posthoctemptestfull$t.ratio >= 36)] <- 36
            
            uniqueconstants <- c('')
            if (factorsinvolvedL > 1) {
              uniqueconstants <- unique(posthoctemptestfull[,2])
            }
            for (cUConstant in 1:length(uniqueconstants)) {
              
              if (factorsinvolvedL > 1) {
                decompconst <- sprintf('%s_when_%s_is_%s', otherfactorsinvolved[1], currentfactor[1], uniqueconstants[cUConstant])
                posthoctemptestsub <- posthoctemptestfull[which(posthoctemptestfull[,2] == uniqueconstants[cUConstant]),]
              } else {
                decompconst <- currentfactor[1]
                posthoctemptestsub <- posthoctemptestfull
              }
              
              subouttable <- data.frame(matrix(NA, nrow=nrow(posthoctemptestsub), ncol=length(colsofinterest)))
              colnames(subouttable) <- colsofinterest
              subouttable$decomp <- decomptext
              subouttable$hold <- decompconst
              
              # hold current factor level constant and subset
              for (currentContrastLine in 1:nrow(posthoctemptestsub)) {
                
                # restrict data to only the constant
                subworkingdatabase <- tempdbs
                if (factorsinvolvedL > 1) {
                  subworkingdatabase <- subworkingdatabase[which(subworkingdatabase[,currentfactor[1]] == posthoctemptestsub[currentContrastLine,2]),]
                } else {
                  otherfactorsinvolved <- currentfactor[1]
                }
                
                contrastlevels <- stringr::str_split(posthoctemptestsub[currentContrastLine,1], ' - ')[[1]]
                contrastdata <- list()
                contrastdata$C1_name <- contrastlevels[1]
                contrastdata$C2_name <- contrastlevels[2]
                
                checkname1 <- contrastlevels[1]
                checkname2 <- contrastlevels[2]
                # if factor levels are numbers rather than strings an error can occur
                checkC <- length(which(subworkingdatabase[,otherfactorsinvolved[1]] == contrastlevels[1])) + length(which(subworkingdatabase[,otherfactorsinvolved[1]] == contrastlevels[2]))
                if (checkC == 0) {
                  #warning('lmerPosthoc(): the factor %s appears to be a numeric column rather than string.\n')
                  checkname1 <- stringr::str_split(contrastlevels[1], factorsinvolved)[[1]][2]
                  checkname2 <- stringr::str_split(contrastlevels[2], factorsinvolved)[[1]][2]
                }
                
                tempvect <- subworkingdatabase[which(subworkingdatabase[,otherfactorsinvolved[1]] == checkname1), dependentvariable[1]]
                tempvectunique <- subworkingdatabase[which(subworkingdatabase[,otherfactorsinvolved[1]] == checkname1), subjectid[1]]
                tempvect <- tempvect[which(!is.na(tempvect))]
                contrastdata$C1_n <- length(tempvect)
                contrastdata$C1_uniquen <- length(unique(tempvectunique))
                contrastdata$C1_mean <- mean(tempvect, na.rm=TRUE)
                contrastdata$C1_sd <- sd(tempvect, na.rm=TRUE)
                if (factorsinvolvedL > 1) {
                  tempvect <- subworkingdatabase[which(subworkingdatabase[,otherfactorsinvolved[1]] == checkname2), dependentvariable[1]]
                  tempvectunique <- subworkingdatabase[which(subworkingdatabase[,otherfactorsinvolved[1]] == checkname2), subjectid[1]]
                } else {
                  tempvect <- subworkingdatabase[which(subworkingdatabase[,currentfactor[1]] == checkname2), dependentvariable[1]]
                  tempvectunique <- subworkingdatabase[which(subworkingdatabase[,currentfactor[1]] == checkname2), subjectid[1]]
                }
                tempvect <- tempvect[which(!is.na(tempvect))]
                contrastdata$C2_n <- length(tempvect)
                contrastdata$C2_uniquen <- length(unique(tempvectunique))
                contrastdata$C2_mean <- mean(tempvect, na.rm=TRUE)
                contrastdata$C2_sd <- sd(tempvect, na.rm=TRUE)
                
                contrastdata$df <- as.numeric(posthoctemptestsub$df[currentContrastLine])
                contrastdata$t <- as.numeric(posthoctemptestsub$t.ratio[currentContrastLine])
                contrastdata$t <- abs(contrastdata$t)
                contrastdata$p <- as.numeric(posthoctemptestsub$p.value[currentContrastLine])
                
                
                ncp <- pkgcond::suppress_conditions(suppressWarnings(MBESS::conf.limits.nct(ncp = as.numeric(contrastdata$t),
                                                                                            df = as.numeric(contrastdata$df),
                                                                                            conf.level = confidenceinterval)))
                
                subouttable$contrast[currentContrastLine] <- currentfactor[1]
                subouttable$df[currentContrastLine] <- contrastdata$df
                subouttable$t.ratio[currentContrastLine] <- contrastdata$t
                subouttable$t.conf.int.lower[currentContrastLine] <- ncp$Lower.Limit
                subouttable$t.conf.int.upper[currentContrastLine] <- ncp$Upper.Limit
                
                subouttable$p.value[currentContrastLine] <- contrastdata$p
                subouttable$p.conf.int.lower[currentContrastLine] <- 2*pt(ncp$Upper.Limit, contrastdata$df, lower=FALSE) 
                subouttable$p.conf.int.upper[currentContrastLine] <- 2*pt(ncp$Lower.Limit, contrastdata$df, lower=FALSE) 
                if (subouttable$p.conf.int.upper[currentContrastLine] > 0.99) {
                  subouttable$p.conf.int.upper[currentContrastLine] <- 0.99
                }
                
                #cat(sprintf('p val comparison: %.2f (model), %.2f (calc)\n', contrastdata$p, (2*pt(abs(contrastdata$t), contrastdata$df, lower=FALSE))))
                
                subouttable$C1name[currentContrastLine] <- contrastdata$C1_name
                subouttable$C1n[currentContrastLine] <- contrastdata$C1_n
                subouttable$C1mean[currentContrastLine] <- contrastdata$C1_mean
                subouttable$C1sd[currentContrastLine] <- contrastdata$C1_sd
                
                subouttable$C2name[currentContrastLine] <- contrastdata$C2_name
                subouttable$C2n[currentContrastLine] <- contrastdata$C2_n
                subouttable$C2mean[currentContrastLine] <- contrastdata$C2_mean
                subouttable$C2sd[currentContrastLine] <- contrastdata$C2_sd
                
                outPvalue <- fuzzyP(contrastdata$p)
                subouttable$significant[currentContrastLine] <- outPvalue$significance
                
                # effect size
                if (otherfactorsinvolved %in% within) {
                  # within subjects
                  
                  # collapse across unnecesary levels - no need given ML approach
                  tempcal <- sprintf("subworkingdatabase <- doBy::summaryBy(%s ~ %s + %s, FUN=c(mean), data=subworkingdatabase, keep.names=TRUE)", dependentvariable[1], subjectid[1], paste(factorsinvolved, collapse=sprintf("+")))
                  suppressWarnings(eval(parse(text=tempcal)))
                  
                  # obtain correlation estimate
                  uniqueids <- unique(subworkingdatabase[,subjectid[1]])
                  cortestdata <- data.frame(matrix(NA, nrow=length(uniqueids), ncol=2))
                  colnames(cortestdata) <- c('C1', 'C2')
                  for (cRCorrtest in 1:length(uniqueids)) {
                    checkindx <- which(subworkingdatabase[,subjectid[1]] == uniqueids[cRCorrtest] &
                                         subworkingdatabase[,otherfactorsinvolved[1]] == contrastdata$C1_name)
                    if (length(checkindx) > 0) {
                      cortestdata[cRCorrtest, 1] <- subworkingdatabase[checkindx, dependentvariable[1]]
                    }
                    checkindx <- which(subworkingdatabase[,subjectid[1]] == uniqueids[cRCorrtest] &
                                         subworkingdatabase[,otherfactorsinvolved[1]] == contrastdata$C2_name)
                    if (length(checkindx) > 0) {
                      cortestdata[cRCorrtest, 2] <- subworkingdatabase[checkindx, dependentvariable[1]]
                    }
                  }
                  
                  # stupid simple hack for insufficient finite observations
                  for (cRCorrtest in 1:4) {
                    if ((length(cortestdata[which(!is.na(cortestdata[,1])),1]) < 5) | (length(cortestdata[which(!is.na(cortestdata[,2])),2]) < 5)) {
                      cortestdata <- rbind(cortestdata, cortestdata, cortestdata)
                    }
                  }
                  correlationtestresult <- tryCatch({
                    correlationtestresult <- tryCatch({
                      correlationtest <- stats::cor.test(cortestdata$C1, cortestdata$C2, alternative='two.sided', method = "pearson", conf.level = confidenceinterval, use = "complete.obs")
                      correlationtestresult <- correlationtest$estimate[[1]]
                      
                      # Percentage bend - robust to outlier
                      sR2 <- WRS2::pbcor(cortestdata$C1, cortestdata$C2, beta=0.2, ci=FALSE, alpha=1-confidenceinterval)
                      correlationtestresult <- sR2$cor
                    }, error = function(e) {
                      correlationtest <- stats::cor.test(cortestdata$C1, cortestdata$C2, alternative='two.sided', method = "pearson", conf.level = confidenceinterval, use = "complete.obs")
                      correlationtestresult <- correlationtest$estimate[[1]]
                    })
                  }, error = function(e) {
                    correlationtestresult <- 0.5 # assumption
                  })
                  subouttable$correlation[currentContrastLine] <- correlationtestresult
                  
                  contrastdata$effectsize <- contrastdata$t * sqrt((2*(1-correlationtestresult))/nrow(cortestdata))
                  contrastdata$effectsize.conf.int.lower <- ncp$Lower.Limit * sqrt((2*(1-correlationtestresult))/nrow(cortestdata))
                  contrastdata$effectsize.conf.int.upper <- ncp$Upper.Limit * sqrt((2*(1-correlationtestresult))/nrow(cortestdata))
                } else {
                  # between subjects
                  contrastdata$effectsize <- contrastdata$t * sqrt((1/contrastdata$C1_uniquen) + (1/contrastdata$C2_uniquen))
                  contrastdata$effectsize.conf.int.lower <- ncp$Lower.Limit * sqrt((1/contrastdata$C1_uniquen) + (1/contrastdata$C2_uniquen))
                  contrastdata$effectsize.conf.int.upper <- ncp$Upper.Limit * sqrt((1/contrastdata$C1_uniquen) + (1/contrastdata$C2_uniquen))
                }
                
                subouttable$effectsize[currentContrastLine] <- contrastdata$effectsize
                subouttable$effectsize.conf.int.lower[currentContrastLine] <- contrastdata$effectsize.conf.int.lower
                subouttable$effectsize.conf.int.upper[currentContrastLine] <- contrastdata$effectsize.conf.int.upper
                
                idtextout <- sprintf('posthoc-%s-%s-vs-%s', factortag, contrastdata$C1_name, contrastdata$C2_name)
                subouttable$idtag[currentContrastLine] <- idtextout
                
              } # currentContrastLine
              
              outtable <- rbind(outtable, subouttable)
              
            } # hold constant
          } # error trap for emmeans
          
          # need to output table to res
          if (nrow(outtable) > 0) {
            textcall <- sprintf("res$Posthoc_%s <- outtable", factortag)
            eval(parse(text=textcall))
          }
        } # no anova required
        
        if (progressbar) {
          # update progress
          #  cat(sprintf('.'))
          pb$tick()
          countticks <- countticks + 1
        }
      }
      if (progressbar) {
        if (countticks < stepcounts) {
          pb$tick(stepcounts)
        }
      }
    } # posthoc limit
  }
  
  return(res)
}