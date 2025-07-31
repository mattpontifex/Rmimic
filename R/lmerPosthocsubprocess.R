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
#' @importFrom stats cor.test formula model.frame complete.cases
#' @importFrom WRS2 pbcor
#' @importFrom data.table data.table as.data.table
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
            subworkingdatabase <- data.table::as.data.table(subworkingdatabase)
            subworkingdatabase <- subworkingdatabase[get(currentfactor[1]) == currentfactorlevelsinvolved[cD]]
            
            currentfactorlevelstring <- currentfactorlevelsinvolved[cD]
            # see if factor is just a number
            if (!grepl("\\D", currentfactorlevelstring)) {
              currentfactorlevelstring <- sprintf('%s%s', currentfactor[1], currentfactorlevelstring)
            }
            decompconst <- sprintf('When %s is %s', currentfactor[1], currentfactorlevelstring)
            
            # validate if random factors can be included
            rf_str <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
            rf_terms <- stringr::str_split(stringr::str_remove_all(rf_str, " "), "~")[[1]][2]
            random_terms <- trimws(unlist(strsplit(rf_terms, "\\+")))
            keep_terms <- character(0)
            for (term in random_terms) {
              # Extract variable names
              vars <- unlist(regmatches(term, gregexpr("\\b\\w+\\b", term)))
              # Only keep variables that are columns in the data
              valid_vars <- intersect(vars, colnames(subworkingdatabase))
              # Vectorized unique count check
              uniq_counts <- sapply(valid_vars, function(v) data.table::uniqueN(subworkingdatabase[[v]]))
              # For all variables in term, must have at least 2 levels
              if (all(uniq_counts > 1)) {
                keep_terms <- c(keep_terms, term)
              } else {
                # Try replacing bad variables with 1 and clean up
                fixed_term <- term
                for (v in valid_vars[uniq_counts <= 1]) {
                  fixed_term <- stringr::str_replace_all(fixed_term, v, "1")
                  fixed_term <- stringr::str_replace_all(fixed_term, ":1|1:", "")
                }
                # If after replacement it's not trivial, keep it
                if (!fixed_term %in% keep_terms && fixed_term != "(1|1)" && fixed_term != "") {
                  keep_terms <- c(keep_terms, fixed_term)
                }
                # Otherwise, discard the term
              }
            }
            randomformula <- paste(keep_terms[keep_terms != ""], collapse = " + ")
            
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
              
            # run the test on the factor of interest using the subsetted data
            smp <- subworkingdatabase
            subfit <- tryCatch({
              model_formula <- as.formula(sprintf('%s ~ %s + %s', dependentvariable[1], fixedformula, randomformula))
              subfit <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(lmerTest::lmer(formula = model_formula, data = smp)))))
            }, error = function(e) {
              cat(sprintf('\nThis error is usually the result of too many random factors specified for the interaction.\nUnable to run decomposition anova with the formula: %s ~ %s + %s\n', dependentvariable[1], fixedformula, randomformula))
              subfit <- NULL
            })
            if (!is.null(subfit)) {
              # get table
              subresults <- invisible(suppressWarnings(suppressMessages(lmerEffects(subfit, dependentvariable=dependentvariable, subjectid = subjectid, within=within, df = df, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, suppresstext=FALSE, smp=smp))))
              # get posthoc
              subresults <- invisible(suppressWarnings(suppressMessages(lmerPosthoc(subresults, between=between, within=within, covariates=covariates, planned=planned, posthoccorrection='none', confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, posthoclimit=posthoclimit-1, calltype='subprocess', progressbar=progressbar))))
              
              # need to output subresults to res
              if (length(names(subresults)) > 0) {
                res[[sprintf("ANOVA_%s", factortag)]] <- subresults
              }
            } # submodel
          } # each factor level involved
        } else {
          # compute all pairwise contrasts
          posthoctemptestfull <- tryCatch({
            emm_formula <- as.formula(paste("~", paste(effectofinterest, collapse = "*")))
            if (factorsinvolvedL > 1) {
              posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fit, emm_formula, adjust = 'none', mode=emmeansdf), by=currentfactor, adjust='none'))))))
            } else {
              posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fit, emm_formula, adjust = 'none', mode=emmeansdf), adjust='none'))))))
            }
            posthoctemptestfull <- as.data.frame(posthoctemptestfull)
            
            # rank deficiencies and when the estimability package identifies an interaction will result in NA being returned
            checkindx <- which(is.na(posthoctemptestfull$estimate))
            if (length(checkindx) > 0) {
              # not all contrasts were returned
              
              # recompute a model with only the effect of interest
              randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
              randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
              fixedformula <- paste(effectofinterest, collapse=sprintf("*"))
              if (factorsinvolvedL > 1) {
                fixedformula <- paste(c(fixedformula, currentfactor), collapse=sprintf("*"))
              }
              
              model_formula <- as.formula(sprintf('%s ~ %s + %s', dependentvariable[1], fixedformula, randomformula[2]))
              fixfit <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(lmerTest::lmer(formula = model_formula, data = smp)))))
              
              # now compute the pairwise contrast
              if (factorsinvolvedL > 1) {
                posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fit, emm_formula, adjust = 'none', mode=emmeansdf), by=currentfactor, adjust='none'))))))
              } else {
                posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fit, emm_formula, adjust = 'none', mode=emmeansdf), adjust='none'))))))
              }
              posthoctemptestfull <- as.data.frame(posthoctemptestfull)
            }
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
                subworkingdatabase <- data.table::as.data.table(subworkingdatabase)
                if (factorsinvolvedL > 1) {
                  subworkingdatabase <- subworkingdatabase[get(currentfactor[1]) == posthoctemptestsub[currentContrastLine, 2]]
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
                checkC <- subworkingdatabase[get(otherfactorsinvolved[1]) %in% contrastlevels, .N]
                if (checkC == 0) {
                  #warning('lmerPosthoc(): the factor %s appears to be a numeric column rather than string.\n')
                  checkname1 <- stringr::str_split(contrastlevels[1], factorsinvolved)[[1]][2]
                  checkname2 <- stringr::str_split(contrastlevels[2], factorsinvolved)[[1]][2]
                }
                
                
                tempdt1 <- subworkingdatabase[get(otherfactorsinvolved[1]) == checkname1]
                tempvect1 <- tempdt1[[dependentvariable[1]]]
                tempvectunique1 <- tempdt1[[subjectid[1]]]
                tempvect1 <- tempvect1[!is.na(tempvect1)]
                
                contrastdata$C1_n <- length(tempvect1)
                contrastdata$C1_uniquen <- data.table::uniqueN(tempvectunique1)
                contrastdata$C1_mean <- mean(tempvect1)
                contrastdata$C1_sd <- sd(tempvect1)
                
                if (factorsinvolvedL > 1) {
                  tempdt2 <- subworkingdatabase[get(otherfactorsinvolved[1]) == checkname2]
                } else {
                  tempdt2 <- subworkingdatabase[get(currentfactor[1]) == checkname2]
                }
                tempvect2 <- tempdt2[[dependentvariable[1]]]
                tempvectunique2 <- tempdt2[[subjectid[1]]]
                tempvect2 <- tempvect2[!is.na(tempvect2)]
                
                contrastdata$C2_n <- length(tempvect2)
                contrastdata$C2_uniquen <- data.table::uniqueN(tempvectunique2)
                contrastdata$C2_mean <- mean(tempvect2)
                contrastdata$C2_sd <- sd(tempvect2)
                
                contrastdata$df <- as.numeric(posthoctemptestsub$df[currentContrastLine])
                contrastdata$t <- as.numeric(posthoctemptestsub$t.ratio[currentContrastLine])
                contrastdata$t <- abs(contrastdata$t)
                contrastdata$p <- as.numeric(posthoctemptestsub$p.value[currentContrastLine])
                
                
                ncp <- invisible(pkgcond::suppress_conditions(suppressMessages(suppressWarnings(MBESS::conf.limits.nct(ncp = as.numeric(contrastdata$t),
                                                                                            df = as.numeric(contrastdata$df),
                                                                                            conf.level = confidenceinterval)))))
                
                subouttable$contrast[currentContrastLine] <- currentfactor[1]
                subouttable$df[currentContrastLine] <- contrastdata$df
                subouttable$t.ratio[currentContrastLine] <- contrastdata$t
                subouttable$t.conf.int.lower[currentContrastLine] <- ncp$Lower.Limit
                subouttable$t.conf.int.upper[currentContrastLine] <- ncp$Upper.Limit
                
                subouttable$p.value[currentContrastLine] <- contrastdata$p
                subouttable$p.conf.int.lower[currentContrastLine] <- 2*pt(ncp$Upper.Limit, contrastdata$df, lower=FALSE) 
                subouttable$p.conf.int.upper[currentContrastLine] <- min(c(2*pt(ncp$Lower.Limit, contrastdata$df, lower=FALSE),
                                                                           0.99)) 

                
                #cat(sprintf('p val comparison: %.2f (model), %.2f (calc)\n', contrastdata$p, (2*pt(abs(contrastdata$t), contrastdata$df, lower=FALSE))))
                
                subouttable$C1name[currentContrastLine] <- contrastdata$C1_name
                subouttable$C1n[currentContrastLine] <- contrastdata$C1_n
                subouttable$C1mean[currentContrastLine] <- contrastdata$C1_mean
                subouttable$C1sd[currentContrastLine] <- contrastdata$C1_sd
                
                subouttable$C2name[currentContrastLine] <- contrastdata$C2_name
                subouttable$C2n[currentContrastLine] <- contrastdata$C2_n
                subouttable$C2mean[currentContrastLine] <- contrastdata$C2_mean
                subouttable$C2sd[currentContrastLine] <- contrastdata$C2_sd
                
                outPvalue <- fuzzyP(contrastdata$p, studywiseAlpha=studywiseAlpha)
                subouttable$significant[currentContrastLine] <- outPvalue$significance
                
                # effect size
                if (otherfactorsinvolved %in% within) {
                  # within subjects
                  
                  # collapse across unnecesary levels
                  subworkingdatabase <- data.table::as.data.table(subworkingdatabase)
                  bycols <- c(subjectid[1], factorsinvolved)
                  subworkingdatabase <- subworkingdatabase[, lapply(.SD, mean, na.rm=TRUE), 
                                                           by = bycols, 
                                                           .SDcols = dependentvariable[1]]
                  
                  # Get data for C1 and C2
                  c1data <- subworkingdatabase[get(otherfactorsinvolved[1]) == contrastdata$C1_name, 
                                               .(C1 = get(dependentvariable[1])), 
                                               by = eval(subjectid[1])]
                  
                  c2data <- subworkingdatabase[get(otherfactorsinvolved[1]) == contrastdata$C2_name, 
                                               .(C2 = get(dependentvariable[1])), 
                                               by = eval(subjectid[1])]
                  cortestdata <- merge(c1data, c2data, by = subjectid[1])
                  cortestdata <- cortestdata[stats::complete.cases(cortestdata),]
                  cortestdata <- as.data.frame(cortestdata)
                  
                  # stupid simple hack for insufficient finite observations
                  for (cRCorrtest in 1:4) {
                    if (nrow(cortestdata) < 5) {
                      cortestdata <- rbind(cortestdata, cortestdata, cortestdata)
                    }
                  }
                  correlationtestresult <- tryCatch({
                    correlationtestresult <- tryCatch({
                      # Percentage bend - robust to outlier
                      sR2 <- invisible(suppressWarnings(suppressMessages(WRS2::pbcor(cortestdata$C1, cortestdata$C2, beta=0.2, ci=FALSE, alpha=1-confidenceinterval))))
                      correlationtestresult <- sR2$cor
                      if (is.na(correlationtestresult)) {
                        # if fails to return a result, run standard correlation
                        correlationtest <- stats::cor.test(cortestdata$C1, cortestdata$C2, alternative='two.sided', method = "pearson", conf.level = confidenceinterval, use = "complete.obs")
                        correlationtestresult <- correlationtest$estimate[[1]]
                      }
                      correlationtestresult <- correlationtestresult
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
            res[[sprintf("Posthoc_%s", factortag)]] <- outtable
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