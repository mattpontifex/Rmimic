#' lmerPosthocsubprocessContrasts
#'
#' @description Subfunction of Rmimic::lmerPosthoc. Performs posthoc test contrasts.

#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, August 19, 2025
#'
#' @importFrom stringr str_split str_remove_all
#' @importFrom pkgcond suppress_conditions
#' @importFrom lmerTest lmer
#' @importFrom emmeans emmeans
#' @importFrom MBESS conf.limits.nct
#' @importFrom stats cor.test formula model.frame complete.cases
#' @importFrom WRS2 pbcor
#' @importFrom data.table data.table as.data.table
#' 
#' @export

lmerPosthocsubprocessContrasts <- function(tempfit, emmeansdf, effectofinterest, currentfactor, factorsinvolved, otherfactorsinvolved, dependentvariable, decomptext, subjectid, confidenceinterval, studywiseAlpha, within, factortag) {

  tempdbs <- stats::model.frame(tempfit)
  
  # establish output
  colsofinterest <- c("decomp", "hold","contrast", "df",
                      "t.ratio", 't.conf.int.lower', 't.conf.int.upper',
                      "p.value", 'p.conf.int.lower', 'p.conf.int.upper',
                      'effectsize', 'effectsize.conf.int.lower', 'effectsize.conf.int.upper', 'significant', 'correlation',
                      'C1name', 'C1n', 'C1mean', 'C1sd',
                      'C2name', 'C2n', 'C2mean', 'C2sd', 'idtag')
  outtable <- data.frame(matrix(NA, nrow=0, ncol=length(colsofinterest)))
  colnames(outtable) <- colsofinterest
  
  # compute all pairwise contrasts
  posthoctemptestfull <- tryCatch({
    emm_formula <- as.formula(paste("~", paste(effectofinterest, collapse = "*")))
    if (length(factorsinvolved) > 1) {
      posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(tempfit, emm_formula, adjust = 'none', mode=emmeansdf), by=currentfactor, adjust='none'))))))
    } else {
      posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(tempfit, emm_formula, adjust = 'none', mode=emmeansdf), adjust='none'))))))
    }
    posthoctemptestfull <- as.data.frame(posthoctemptestfull)
    
    # rank deficiencies and when the estimability package identifies an interaction will result in NA being returned
    checkindx <- which(is.na(posthoctemptestfull$estimate))
    if (length(checkindx) > 0) {
      # not all contrasts were returned
      
      # recompute a model with only the effect of interest
      randomformula <- Reduce(paste, deparse(stats::formula(tempfit, random.only = TRUE)))
      randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
      fixedformula <- paste(effectofinterest, collapse=sprintf("*"))
      if (length(factorsinvolved) > 1) {
        fixedformula <- paste(c(fixedformula, currentfactor), collapse=sprintf("*"))
      }
      
      model_formula <- as.formula(sprintf('%s ~ %s + %s', dependentvariable[1], fixedformula, randomformula[2]))
      fixfit <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(lmerTest::lmer(formula = model_formula, data = tempdbs)))))
      
      # now compute the pairwise contrast
      if (length(factorsinvolved) > 1) {
        posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fixfit, emm_formula, adjust = 'none', mode=emmeansdf), by=currentfactor, adjust='none'))))))
      } else {
        posthoctemptestfull <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(summary(pairs(emmeans::emmeans(fixfit, emm_formula, adjust = 'none', mode=emmeansdf), adjust='none'))))))
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
    if (length(factorsinvolved) > 1) {
      uniqueconstants <- unique(posthoctemptestfull[,2])
    }
    for (cUConstant in 1:length(uniqueconstants)) {
      
      if (length(factorsinvolved) > 1) {
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
        if (length(factorsinvolved) > 1) {
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
        
        if (length(factorsinvolved) > 1) {
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
  
  return(outtable)
}
  