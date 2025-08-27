#' lmerEffects
#'
#' @description Compute SPSS style univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. 
#'
#' @param fit modelfit from the lmer function
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label. If left NULL assumes all factors are between subjects
#' @param within text list specifying any variables that should be considered as within subjects
#' @param df Parameter to indicate what degrees of freedom approximation should be used. Default is Kenward-Roger. Other options are Shattertwaite or Traditional.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param suppresstext Boolean parameter to avoid obtaining html formatted text results.
#'
#' @return
#' \item{stats}{ANOVA summary table.}
#' \item{randomstats}{ANOVA summary table for random effects.}
#' \item{rsquared}{summary table for r squared portions.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom stringr str_split str_remove_all str_replace_all
#' @importFrom utils packageDate
#' @importFrom pkgcond suppress_conditions
#' @importFrom stats anova residuals model.frame formula
#' @importFrom psychometric CI.Rsq
#' @importFrom lme4 isSingular
#' @importFrom lmerTest ranova
#' @importFrom MuMIn r.squaredGLMM
#' 
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'
#' @export

lmerEffects <- function(fit, dependentvariable=NULL, subjectid=NULL, within=NULL, df = NULL, confidenceinterval=0.95, studywiseAlpha=0.05, suppresstext=FALSE, smp=NULL, verbose=FALSE, ...) {
  
  # debug
  debug <- FALSE
  if (debug == TRUE) {
    data = Rmimic::alertness
    data <- data[which(data$Condition == 'Condition2'),]
    subjectid = "PartID"
    df = NULL
    confidenceinterval=0.95
    studywiseAlpha=0.05
    fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = data)
    
    data <- Rmimic::alertness
    data <- data[which(data$Condition == 'Condition2'),]
    #data$NewGroup <- NA
    #data$NewGroup[which(data$Group == 'Group1')] <- 1
    #data$NewGroup[which(data$Group == 'Group2')] <- 2
    #data$Group <- as.numeric(data$NewGroup)
    fit <- lmerTest::lmer(Alertness ~ Group*Time + Sex + (1 | PartID), data = data)
    
    subjectid = "PartID"
    dependentvariable = NULL
    within=NULL
    randomslope=NULL
    df = NULL
  }
  
  options(contrasts = c("contr.sum", "contr.poly"))
  oldw <- getOption("warn")
  options(warn = -1) #
  
  res <- list()
  
  if (is.null(studywiseAlpha)) {
    studywiseAlpha <- 0.05
  }
  if (is.null(confidenceinterval)) {
    confidenceinterval <- 0.95
  }
  
  randomtext <- NULL
  fixedterms <- NULL
  
  tempdbs <- stats::model.frame(fit)
  formulaextraction <- tryCatch({
    fullformula <- Reduce(paste, deparse(stats::formula(fit)))
    fixedformula <- Reduce(paste, deparse(stats::formula(fit, fixed.only = TRUE)))
    randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
    
    fullformula <- stringr::str_split(stringr::str_remove_all(fullformula, ' '), '~')[[1]]
    fixedformula <- stringr::str_split(stringr::str_remove_all(fixedformula, ' '), '~')[[1]]
    randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
    formulaextraction <- 1
  }, error = function(e) {
    formulaextraction <- NULL
  })
  if (!is.null(formulaextraction)) {
    if (is.null(dependentvariable)) {
      # if we dont know what the DV is
      dependentvariable <- tryCatch({
        dependentvariable <- Reduce(intersect, list(fullformula,fixedformula,randomformula))
      }, error = function(e) {
        dependentvariable <- NULL
      })
    }
    if (!is.null(dependentvariable)) {
      fixedformula <- fixedformula[which(fixedformula != dependentvariable)]
      randomformula <- randomformula[which(randomformula != dependentvariable)]
      fixedterms <- stringr::str_split(fixedformula, "\\+|\\*|\\/|[|]")[[1]]
      randomformula <- stringr::str_split(randomformula, "\\+|\\*|\\/")[[1]]
      randomterms <- colnames(tempdbs)
      randomterms <- randomterms[which(randomterms != dependentvariable)]
      randomterms <- randomterms[which(!(randomterms %in% fixedterms))]
      
      randomtext <- c()
      for (cRT in 1:length(randomterms)) {
        for (cRF in 1:length(randomformula)) {
          if (sprintf('(1|%s)', randomterms[cRT]) == randomformula[cRF]) {
            randomtext <- c(randomtext, sprintf('random intercept for %s', randomterms[cRT]))
            randomformula[cRT] <- ''
          }
        }
      }
      randomformula <- randomformula[which(randomformula != '')]
      if (length(randomformula) > 0) {
        for (cRF in 1:length(randomformula)) {
          tempstr <- stringr::str_replace_all(stringr::str_split(randomformula[cRF], "[|]")[[1]], '[(]|[)]', "")
          if (length(tempstr) == 2)
            randomtext <- c(randomtext, sprintf('random intercept and random slope for %s within %s', tempstr[1], tempstr[2]))
        }
      }
    }
  }
  
  outstring <- tryCatch({
    outstring <- ''
    if (!is.null(dependentvariable)) {
      outstring <- sprintf('Analysis of %s were conducted using a', dependentvariable)
      
      fixedformula <- stringr::str_replace_all(fixedformula, '[+]', ' + ') # add space around addition
      fixedformula <- stringr::str_replace_all(fixedformula, '[-]', ' - ') # add space around subtraction
      fixedformula <- stringr::str_replace_all(fixedformula, '[*]', ' \u00D7 ') # swap out multiplication
      if (!is.null(fixedterms)) {
        for (cB in 1:length(fixedterms)) {
          factorname <- tolower(fixedterms[cB])
          factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
          factorsinvolved <- tolower(unique(unlist(as.character(tempdbs[,fixedterms[cB]]))))
          tempstring <- sprintf("%d (%s: %s)", length(factorsinvolved), factorname, paste(factorsinvolved, collapse = ", "))
          fixedformula <- stringr::str_replace_all(fixedformula, sprintf('%s', fixedterms[cB]), tempstring)
        }
      }
      outstring <- sprintf('%s %s', outstring, fixedformula)
      
      outstring <- sprintf('%s univariate', outstring)
      if (!is.null(within)) {
        outstring <- sprintf('%s repeated measures', outstring)
      }
      outstring <- sprintf('%s multi-level model', outstring) 
      
      if (!is.null(randomtext)) {
        outstring <- sprintf("%s including the", outstring)
        for (cB in 1:length(randomtext)) {
          outstring <- sprintf('%s %s', outstring, randomtext[cB])
          if (length(randomtext) > 1) {
            # more than one random term
            if (cB < length(randomtext)) {
              outstring <- sprintf("%s,", outstring)
            }
            if (cB == (length(randomtext)-1)) {
              outstring <- sprintf("%s and", outstring)
            }
          }
        }
      }
      outstring <- sprintf("%s.", outstring)
    
    }
    
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s Analyses were conducted using the lme4 (Bates et al., %s), lmerTest (Kuznetsova et al., %s), and Rmimic (Pontifex, 2020) packages in %s (%s).', outstring, 
                         strsplit(as.character(utils::packageDate("lme4")),"-")[[1]][1],
                         strsplit(as.character(utils::packageDate("lmerTest")),"-")[[1]][1],
                         rvers, gsub("^\\((\\d{4}).*", "\\1", unlist(strsplit(R.version.string, " "))[4]))
    
    
    
    
    
  }, error = function(e) {
    outstring <- NULL
  })
  
  outmessage <- tryCatch({
    if (!is.null(fit@optinfo$conv$lme4$messages)) {
      outmessage <- sprintf('\n\nMessage from lmer: %s.', paste(fit@optinfo$conv$lme4$messages, collapse=". "))
    } else {
      outmessage <- ""
    }
  }, error = function(e) {
    outmessage <- NULL
  })
  if (!is.null(outmessage)) {
    if (!is.null(outstring)) {
      outstring <- sprintf('%s %s', outstring, outmessage)
    } else {
      outstring <- sprintf('%s', outmessage)
    }
  }
  
  starttime <- Sys.time()
  # see if the model was fit with numeric factors
  boolcheckfactors <- rep_len(1, length(colnames(tempdbs)))
  for (cR in 1:length(boolcheckfactors)) {
    if (is.numeric(tempdbs[,cR])) {
      if (!(colnames(tempdbs)[cR] == dependentvariable)) {
        # numeric but not a dependent variable
        warning(sprintf('lmerEffects() the specified model included a numeric factor (%s). Numeric factor levels may result in errors. Attempting to resolve.\n', colnames(tempdbs)[cR]), immediate=TRUE)
        tempdbs[,cR] <- paste0(colnames(tempdbs)[cR], tempdbs[,cR], sep="")
        boolcheckfactors[cR] <- 0
      }
    } 
  }
  if (sum(boolcheckfactors) < length(colnames(tempdbs))) {
    # refit model
    fit <- tryCatch({
      model_formula <- as.formula(Reduce(paste, deparse(stats::formula(fit))))
      fit <- invisible(pkgcond::suppress_conditions(lmerTest::lmer(formula = model_formula, data = tempdbs)))
    }, error = function(e) {
      cat(sprintf('\n\nlmerEffects() Warning: the specified model includes numeric factor(s). Numeric factor levels may result in errors.\n'))
      model_formula <- as.formula(Reduce(paste, deparse(stats::formula(fit))))
      fit <- lmerTest::lmer(formula = model_formula, data = tempdbs)
    })
  }
  if (verbose) {
    cat(sprintf('  lmerEffects(): time to check model - %.2f sec\n', Sys.time() - starttime))
  }
  
  starttime <- Sys.time()
  res$fit <- fit
  if (!is.null(outstring)) {
    if (!(outstring == '')) {
      res$messageout <- outstring
    }
  }
  
  if (!is.null(dependentvariable)) {
    res$dependentvariable <- dependentvariable
  }
  if (!is.null(fixedterms)) {
    res$fixedterms <- fixedterms
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
  res$df <- df
  res$confidenceinterval <- confidenceinterval
  res$studywiseAlpha <- studywiseAlpha

  # check how many unique participants there are
  if (!is.null(subjectid)) {
    res$subjectid <- subjectid
    numparticipants <- tryCatch(
      {
        tempdbs <- data.table::as.data.table(stats::model.frame(fit))
        data.table::uniqueN(tempdbs[[subjectid]])
      },
      error = function(e) NULL
    )
  } else {
    numparticipants <- tryCatch(
      nrow(data.table::as.data.table(stats::model.frame(fit))),
      error = function(e) NULL
    )
  }
  res$numparticipants <- numparticipants
  
  # check how many factors there are
  numfactors <- tryCatch({
    tempdbs <- stats::model.frame(fit)
    numfactors <- length(colnames(tempdbs)) - 2 # remove subject, remove dependent variable
    if (numparticipants == nrow(tempdbs)) {
      numfactors <- numfactors + 1 # no subject level
    }
    numfactors <- as.numeric(numfactors)
  }, error = function(e) {
    numfactors <- NULL
  })
  res$numfactors <- numfactors
  
  # see if the model may be over-fit
  res$booloverfitwarning <- FALSE
  singulartest <- invisible(suppressMessages(suppressWarnings(lme4::isSingular(fit))))
  if (singulartest == TRUE) {
    res$booloverfitwarning <- TRUE
    #warning("Warning: Model may be over-fit.")
  }
  if (verbose) {
    cat(sprintf('  lmerEffects(): time to determine parameters - %.2f sec\n', Sys.time() - starttime))
  }
  
  
  #Time difference of 0.989 secs
  starttime <- Sys.time()
  # Compute ANOVA table for model
  # note KR is computationally intensive and slow for large datasets
  # could use 
  # as <- afex::mixed(as.formula(Reduce(paste, deparse(stats::formula(fit)))),
  #           data = stats::model.frame(fit), method = "KR")   # method = "S" for Satterthwaite
  # however it does not provide sum of square information which means no parial eta squared calculation
  # solution is to use Satterthwaite for large datasets
  as <- tryCatch({
    if (toupper(df) == toupper("Shattertwaite")) {
      as <- invisible(suppressWarnings(suppressMessages(data.frame(pkgcond::suppress_conditions(stats::anova(fit, type = 3))))))
    } else {
      as <- invisible(suppressWarnings(suppressMessages(data.frame(pkgcond::suppress_conditions(stats::anova(fit, type = 3, ddf = "Kenward-Roger"))))))
    }
  }, error = function(e) {
    as <- NULL
  })
  if (verbose) {
    cat(sprintf('  lmerEffects(): time to compute anova - %.2f sec\n', Sys.time() - starttime))
  }
  
  # Supplement ANOVA table for model
  if (!is.null(as)) {
    dataframeout <- data.frame(matrix(NA,nrow=nrow(as),ncol=12))
    colnames(dataframeout) <- c("Effect", "DFn", "DFd", "SSn", "SSd", "SSe", "F.value", "p.value", "partialetasquared", "fsquared","fsquared.ci.lower", "fsquared.ci.upper")
    
    dataframeout['Effect'] <- rownames(as)
    dataframeout['DFn'] <- as.integer(as$NumDF)
    dataframeout['DFd'] <- as.integer(as$DenDF)
    dataframeout['SSn'] <- as$Sum.Sq
    dataframeout['SSd'] <- as$Mean.Sq
    tempvect <- tryCatch({
      tempvect <- sum(dataframeout$SSd, na.rm=TRUE) + mean(stats::residuals(fit)^2, na.rm=TRUE)
    }, error = function(e) {
      tempvect <- NULL
      stop('Unable to compute the mean(stats::residuals(fit)^2) for this model. Please verify that the model is correct.')
    })
    dataframeout['SSe'] <- tempvect
    dataframeout['F.value'] <- as$F.value
    dataframeout['F.value.ci.lower'] <- NA
    dataframeout['F.value.ci.upper'] <- NA
    #stats::pf(F,DFn,DFd, lower.tail=FALSE)
    dataframeout['p.value'] <- as$Pr..F.
    dataframeout$significance <- NA
    for (cT in 1:nrow(dataframeout)) {
      dataframeout$significance[cT] <- Rmimic::fuzzyP(dataframeout$p.value[cT], studywiseAlpha)$significance
    }
    dataframeout['p.value.ci.lower'] <- NA
    dataframeout['p.value.ci.upper'] <- NA
    dataframeout$partialetasquared <- dataframeout$SSn / (dataframeout$SSn + dataframeout$SSe)
    
    #Calculate f^2 = partialeta^2 / ( 1 - partialeta^2 )
    dataframeout$fsquared <- dataframeout$partialetasquared/(1-dataframeout$partialetasquared)
    if (!is.null(numparticipants)) {
      for (cT in 1:nrow(dataframeout)) {
        if (!is.na(dataframeout$partialetasquared[cT])) {
          tempvect <- tryCatch({
            temp <- invisible(pkgcond::suppress_conditions(suppressWarnings(psychometric::CI.Rsq(dataframeout$partialetasquared[cT], numparticipants, length(unlist(strsplit(as.character(dataframeout$Effect[cT]),"[:]"))), level = confidenceinterval))))
          }, error = function(e) {
            tempvect <- NULL
          })
          if (!is.null(tempvect)) {
            if (temp$LCL[1]<0) {
              dataframeout$fsquared.ci.lower[cT] <- 0
            } else {
              dataframeout$fsquared.ci.lower[cT] <- (temp$LCL[1]/(1-temp$LCL[1]))
            }
            if (temp$UCL[1] > 1) {
              dataframeout$fsquared.ci.upper[cT] <- Inf
            } else {
              dataframeout$fsquared.ci.upper[cT] <- (temp$UCL[1]/(1-temp$UCL[1]))
            }
          }
        }
      }
    }
    res$stats <- dataframeout
  }
  
  #Time difference of 0.1391 secs
  starttime <- Sys.time()
  #https://www.rensvandeschoot.com/tutorials/lme4/
  # lmerTest::ranova tries to access data from the parent environment rather than data through the call
  # bootstrap function can pass data through smp 
  as <- tryCatch({
    as <- invisible(data.frame(lmerTest::ranova(fit)))
  }, error = function(e) {
    # oldschool approach
    tempdbs <- stats::model.frame(fit) # pull data out
    fit <- invisible(update(fit, data=tempdbs, evaluate = TRUE)) # tell model to update itself with the old data
    as <- invisible(data.frame(lmerTest::ranova(fit))) # magically this function works again
  })
  if (verbose) {
    cat(sprintf('  lmerEffects(): time to compute ranova - %.2f sec\n', Sys.time() - starttime))
  }
  
  if (!is.null(as)) {
    as <- as[2:nrow(as),]
    dataframeout <- data.frame(matrix(NA,nrow=nrow(as),ncol=5))
    colnames(dataframeout) <- c("Effect", "DF", "LogLikelihood", "LRT", "p.value")
    dataframeout['Effect'] <- rownames(as)
    dataframeout['DF'] <- as$Df
    dataframeout['LogLikelihood'] <- as$logLik
    dataframeout['LRT'] <- as$LRT
    dataframeout['p.value'] <- as$Pr..Chisq.
    dataframeout$significance <- NA
    for (cT in 1:nrow(dataframeout)) {
      dataframeout$significance[cT] <- Rmimic::fuzzyP(dataframeout$p.value[cT], studywiseAlpha)$significance
    }
    res$randomstats <- dataframeout
  }
  
  # Compute R^2
  #Nakagawa, S., & Schielzeth, H. (2012). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4, 133-142. doi: 10.1111/j.2041-210x.2012.00261.x
  #Nakagawa, S., Johnson, P. C., & Schielzeth, H. (2017). The coefficient of determination R^2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface, 14(134), 20170213. doi:10.1098/rsif.2017.0213
  dataframeout <- data.frame(matrix(NA,nrow=3,ncol=4))
  colnames(dataframeout) <- c('portion', 'effects', 'ci.lower', 'ci.upper')
  dataframeout$portion <- c('Fixed', 'Random', 'Model')  
  starttime <- Sys.time()
  as <- tryCatch({
    as <- invisible(suppressWarnings(MuMIn::r.squaredGLMM(fit)))
  }, error = function(e) {
    as <- NULL
    stop('Unable to compute the MuMIn::r.squaredGLMM(fit) for this model. Please verify that the model is correct.')
  })
  if (verbose) {
    cat(sprintf('  lmerEffects(): time to compute rsquared - %.2f sec\n', Sys.time() - starttime))
  }
  if (!is.null(as)) {
    as[3] <- as[2] # model
    #random = model - fixed
    as[2] <- as[3] - as[1]
    ciLower <- rep_len(NA, length(as))
    ciUpper <- rep_len(NA, length(as))
    for (cT in 1:length(as)) {
      temp <- tryCatch({
        temp <- invisible(pkgcond::suppress_conditions(suppressWarnings(psychometric::CI.Rsq(as[cT], numparticipants, numfactors, level=confidenceinterval))))
      }, error = function(e) {
        temp <- NULL
      })
      if (!is.null(temp)) {
        if (temp$LCL[1]<0) {
          temp$LCL[1] <- 0
        } 
        ciLower[cT] <- temp$LCL[1]
        if (temp$UCL[1] > 1) {
          temp$UCL[1] <- Inf
        }
        ciUpper[cT] <- temp$UCL[1]
      }
      dataframeout[cT,2] <- as[cT]
      dataframeout[cT,3] <- ciLower[cT]
      dataframeout[cT,4] <- ciUpper[cT]
    }
    res$rsquared <- dataframeout
  }
  
   starttime <- Sys.time()
  if (!suppresstext) {
    # obtain text outputs
    res <- lmerEffects2text(res)
  }
   if (verbose) {
     cat(sprintf('  lmerEffects(): time to process text- %.2f sec\n', Sys.time() - starttime))
   }
 
  options(warn = oldw) # turn warnings back to original settings
  return(res)
}
  