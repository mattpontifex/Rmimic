#' lmerModelstats
#'
#' @description Obtain standard linear model reporting statistics for a lmer test
#'
#' @param fit lmer object.
#' @param altfit lmer object to use as comparison in hierarchical regression. If NULL then standard contrast against an empty model with just random effects is performed.
#' @param df text parameter to indicate what degrees of freedom to use. Options are Kenward-Roger (default) or Shattertwaite
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#'
#' @return A list of standard reporting statistics
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, June 3, 2025
#'
#' @importFrom stats formula update anova
#' @importFrom stringr str_split str_remove_all
#' @importFrom pkgcond suppress_conditions
#' @importFrom car vif
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom effectsize F_to_eta2
#' @importFrom lme4 fixef
#' @importFrom parameters standardize_parameters
#' 
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerModelstats(fit)
#'     
#'     altfit <- lmerTest::lmer(Alertness ~ Time + (1 | PartID), data = Rmimic::alertness)
#'     fit <- lmerTest::lmer(Alertness ~ Group + Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerModelstats(fit, altfit)
#' 
#' @export

lmerModelstats <- function(fit, altfit = NULL, df = NULL, confidenceinterval=0.95, studywiseAlpha=0.05) {
  
  res <- list()
  
  res$altfit <- altfit
  res$fit <- fit
  res$confidenceinterval <- confidenceinterval
  res$studywiseAlpha <- studywiseAlpha
  
  # determine what degrees of freedom to use
  if (!is.null(df)) {
    if (toupper(df) == toupper("Kenward-Roger")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    }
  } else {
    df = "Kenward-Roger"
  }
  res$dfparam <- df
  
  # create comparison model only containing random effects
  randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
  randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
  subfit <- stats::update(fit, formula(paste("~", randomformula[2])))
  res$dv <- randomformula[1]
  res$fullformula <- Reduce(paste, deparse(stats::formula(fit)))
    
  # Compute ANOVA table for model
  as <- tryCatch({
    if (toupper(df) == toupper("Shattertwaite")) {
      as <- data.frame(pkgcond::suppress_conditions(stats::anova(subfit, fit, type = 3)))
    } else {
      as <- data.frame(pkgcond::suppress_conditions(stats::anova(subfit, fit, type = 3, ddf = "Kenward-Roger")))
    }
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    res$AIC <- as$AIC[2]
    res$df <- as$Df[2]
    res$Chisq <- as$Chisq[2]
    res$F.value <- res$Chisq / res$df
    res$deviance <- as$deviance[2]
    res$p.value <- as$Pr..Chisq.[2]
    pullvalue <- fuzzyP(res$p.value, studywiseAlpha=studywiseAlpha, html=TRUE)
    if (pullvalue$significance) {
      res$significant <- TRUE
    } else {
      res$significant <- FALSE
    }
  }
  
  as <- tryCatch({
    as <- max(car::vif(fit), na.rm=TRUE)
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    res$VIF <- as
  }
  
  # Compute Rsquared
  as <- tryCatch({
    as <- suppressWarnings(MuMIn::r.squaredGLMM(fit))
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    as <- tryCatch({
      if (nrow(as) == 2) {
        as <- as.data.frame(as)
        as <- unlist(as[2,])
      }
      as <- as
    }, error = function(e) {
      as <- as
    })
  }
  if (!is.null(as)) {
    as[3] <- as[2] # model
    #random = model - fixed
    as[2] <- as[3] - as[1]
    res$fixed.r.squared <- as[1]
    res$random.r.squared <- as[2]
    res$model.r.squared <- as[3]
  }
  
  # Compute fsquared
  as <- tryCatch({
    partialetasquared <- effectsize::F_to_eta2(res$F.value, res$df, 1000, ci = confidenceinterval, alternative = "two.sided")
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    res$fsquared <- partialetasquared$Eta2_partial / (1-partialetasquared$Eta2_partial)
    res$fsquared.ci.lower <- partialetasquared$CI_low / (1-partialetasquared$CI_low)
    res$fsquared.ci.upper <- partialetasquared$CI_high / (1-partialetasquared$CI_high)
  }
  
  # Obtain Coefficients
  nam <- names(lme4::fixef(fit))
  varsofinterest <- c('Variable','B','B.lower.conf.int','B.upper.conf.int','SE','Beta','Beta.lower.conf.int','Beta.upper.conf.int','t','z','df','p.value', 'significant','OddsRatio','OddsRatio.lower.conf.int','OddsRatio.upper.conf.int')
  res$coefficients <- data.frame(matrix(data=NA,nrow=length(nam),ncol=length(varsofinterest)))
  names(res$coefficients) <- varsofinterest
  res$coefficients$Variable <- nam
  msc <- as.data.frame(summary(fit)$coefficients)
  if ('Estimate' %in% names(msc)) {res$coefficients$B <- msc$Estimate}
  if ("Std. Error" %in% names(msc)) {res$coefficients$SE <- msc$`Std. Error`}
  if ("df" %in% names(msc)) {res$coefficients$df <- msc$df}
  if ("t value" %in% names(msc)) {res$coefficients$t <- msc$`t value`}
  if ("z value" %in% names(msc)) {res$coefficients$z <- msc$`z value`}
  if ('Pr(>|t|)' %in% names(msc)) {res$coefficients$p.value <- msc$`Pr(>|t|)`}
  if ('Pr(>|z|)' %in% names(msc)) {res$coefficients$p.value <- msc$`Pr(>|z|)`}
  res$coefficients$significant <- FALSE
  for (cR in 1:nrow(res$coefficients)) {
    pullvalue <- Rmimic::fuzzyP(res$coefficients$p.value[cR], studywiseAlpha=studywiseAlpha, html=TRUE)
    if (pullvalue$significance) {
      res$coefficients$significant[cR] <- TRUE
    }
  }
  
  # obtain confidence intervals
  as <- tryCatch({
    as <- as.data.frame(suppressMessages(suppressWarnings(confint(fit, level=confidenceinterval))))
    as <- as[which(rownames(as) %in% nam),]
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    res$coefficients$B.lower.conf.int <- as$`2.5 %`
    res$coefficients$B.upper.conf.int <- as$`97.5 %`
  }
  
  # obtain beta
  as <- tryCatch({
    as <- as.data.frame(parameters::standardize_parameters(fit, ci=confidenceinterval))
  }, error = function(e) {
    as <- NULL
  })
  if (!is.null(as)) {
    res$coefficients$Beta <- as$Std_Coefficient
    res$coefficients$Beta.lower.conf.int <- as$CI_low
    res$coefficients$Beta.upper.conf.int <- as$CI_high
    
    # see if we should be reporting odds ratios instead
    tempvect <- family(fit)
    if (tempvect$family == 'binomial') {
      res$coefficients$OddsRatio <- exp(as$Std_Coefficient)
      res$coefficients$OddsRatio.lower.conf.int <- exp(as$CI_low)
      res$coefficients$OddsRatio.upper.conf.int <- exp(as$CI_high)
    }
  }
  
  
  if (!is.null(altfit)) {
    res$change <- list()
    
    res$altformula <- Reduce(paste, deparse(stats::formula(altfit)))
    as <- tryCatch({
      if (toupper(df) == toupper("Shattertwaite")) {
        as <- data.frame(pkgcond::suppress_conditions(stats::anova(altfit, fit, type = 3)))
      } else {
        as <- data.frame(pkgcond::suppress_conditions(stats::anova(altfit, fit, type = 3, ddf = "Kenward-Roger")))
      }
    }, error = function(e) {
      as <- NULL
    })
    if (!is.null(as)) {
      res$change$AIC <- as$AIC[2]
      res$change$df <- as$Df[2]
      res$change$Chisq <- as$Chisq[2]
      res$change$F.value <- res$Chisq / res$df
      res$change$deviance <- as$deviance[2]
      res$change$p.value <- as$Pr..Chisq.[2]
      pullvalue <- Rmimic::fuzzyP(res$change$p.value, studywiseAlpha=studywiseAlpha, html=TRUE)
      if (pullvalue$significance) {
        res$change$significant <- TRUE
      } else {
        res$change$significant <- FALSE
      }
    }
    
    as <- tryCatch({
      as <- suppressWarnings(MuMIn::r.squaredGLMM(altfit))
    }, error = function(e) {
      as <- NULL
    })
    if (!is.null(as)) {
      as <- tryCatch({
        if (nrow(as) == 2) {
          as <- as.data.frame(as)
          as <- unlist(as[2,])
        }
        as <- as
      }, error = function(e) {
        as <- as
      })
    }
    if (!is.null(as)) {
      as[3] <- as[2] # model
      #random = model - fixed
      as[2] <- as[3] - as[1]
      res$change$fixed.r.squared <- res$fixed.r.squared - as[1]
      res$change$random.r.squared <- res$random.r.squared - as[2]
      res$change$model.r.squared <- res$model.r.squared - as[3]
    }
    
    # Compute fsquared
    as <- tryCatch({
      partialetasquared <- effectsize::F_to_eta2(res$change$F.value, res$change$df, 1000, ci = confidenceinterval, alternative = "two.sided")
    }, error = function(e) {
      as <- NULL
    })
    if (!is.null(as)) {
      res$change$fsquared <- partialetasquared$Eta2_partial / (1-partialetasquared$Eta2_partial)
      res$change$fsquared.ci.lower <- partialetasquared$CI_low / (1-partialetasquared$CI_low)
      res$change$fsquared.ci.upper <- partialetasquared$CI_high / (1-partialetasquared$CI_high)
    }
  }
  
  # generate text outputs
  res <- lmerModelstats2text(res, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha)
    
  # add some process information
  outstring <- ''
  rvers <- unlist(strsplit(R.version.string, " "))
  rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
  outstring <- sprintf('%s Regression analyses were conducted using the lme4 (Bates et al., %s), lmerTest (Kuznetsova et al., %s), and Rmimic (Pontifex, %s) packages in %s.', outstring, 
                       strsplit(as.character(utils::packageDate("lme4")),"-")[[1]][1],
                       strsplit(as.character(utils::packageDate("lmerTest")),"-")[[1]][1],
                       strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1], rvers)
  
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
  
  res$messageout <- outstring

  return(res)
}
