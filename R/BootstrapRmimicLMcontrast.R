#' BootstrapRmimicLMcontrast
#'
#' @description Compute results for bootstrapped regression analysis with effect size and confidence intervals. This function takes stats::lm fits for a base model and the model of interest and calculates statistics for the model of interest relative to the base model. For logistic regression, pseudo r2 values are reported using Tjurs 2009 approach.
#'
#' @param fit fit from a basic linear model
#' @param altfit fit from a linear model of interest
#' @param repetitions number of bootstrap repetitions. Default is 999.
#' @param resample_min minimum number of samples to include when resampling.
#' @param resample_max minimum number of samples to include when resampling.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#'
#' @return
#' \item{stats}{ANOVA summary table.}
#' \item{modelcheck}{chisquared summary table.}
#' \item{changestats}{ANOVA summary table of change statistics.}
#' \item{coefficients}{A summary table of regression coefficients.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, March 25, 2025
#' 
#' @importFrom  dplyr slice_sample
#' @importFrom stats runif
#' 
#'
#' @examples
#'
#'     # Compute Model from mtcars dataset
#'     basefit <- lm(mpg ~ am + wt, data = mtcars)
#'     fit <- lm(mpg ~ am + wt + qsec, data = mtcars)
#'     result <- BootstrapRmimicLMcontrast(basefit, fit,
#'          repetitions=99, resample_min=32, resample_max=96)
#'
#' @export

BootstrapRmimicLMcontrast <- function(basefit, fit, repetitions=999, resample_min=FALSE, resample_max=FALSE, studywiseAlpha=0.05, confidenceinterval=0.95) {
  #repetitions <- 999
  #resample_min <- 32
  #resample_max <- resample_min*3
  #basefit <- stats::lm(mpg ~ 1, data = mtcars)
  #fit <- stats::lm(mpg ~ hp + wt, data = mtcars)
  #result <- BootstrapRmimicLMcontrast(basefit, fit, repetitions=999, resample_min=resample_min, resample_max=resample_max)
  
  if (!resample_min) {
    resample_min <- nrow(model.frame(basefit))
  }
  if (!resample_max) {
    resample_max <- nrow(model.frame(basefit))
  }  
  
  # run model on full dataset
  mainregresult <- Rmimic::RmimicLMcontrast(basefit, fit, studywiseAlpha=studywiseAlpha, confidenceinterval=confidenceinterval, verbose=FALSE)
  
  cat(sprintf('BootstrapRmimicLMcontrast(): Beginning the requested %d model reptitions...\n', repetitions))
  
  # loop through analysis
  for (cN in 1:repetitions) {
    # Determine how large a random sample of the original database
    numberofsamples <- floor(stats::runif(1, min=resample_min, max=resample_max))
    
    # populates dataset subsampled from the original sample with replacement
    smp <- dplyr::slice_sample(model.frame(fit), n=numberofsamples, replace=TRUE)
    
    # run model on this dataset
    basefitnew <- stats::lm(formula(basefit), data=smp)
    fitnew <- stats::lm(formula(fit), data=smp)
    
    # compare models
    regresult <- Rmimic::RmimicLMcontrast(basefitnew, fitnew, studywiseAlpha=studywiseAlpha, confidenceinterval=confidenceinterval, verbose=FALSE)
    
    regresult$stats$Model[which(regresult$stats$Model == "formula(basefit)")] <- Reduce(paste, deparse(formula(basefit)))
    regresult$modelcheck$Model[which(regresult$modelcheck$Model == "formula(basefit)")] <- Reduce(paste, deparse(formula(basefit)))
    regresult$changestats$Model[which(regresult$changestats$Model == "formula(basefit)")] <- Reduce(paste, deparse(formula(basefit)))
    regresult$coefficients$Model[which(regresult$coefficients$Model == "formula(basefit)")] <- Reduce(paste, deparse(formula(basefit)))
    
    regresult$stats$Model[which(regresult$stats$Model == "formula(fit)")] <- Reduce(paste, deparse(formula(fit)))
    regresult$modelcheck$Model[which(regresult$modelcheck$Model == "formula(fit)")] <- Reduce(paste, deparse(formula(fit)))
    regresult$changestats$Model[which(regresult$changestats$Model == "formula(fit)")] <- Reduce(paste, deparse(formula(fit)))
    regresult$coefficients$Model[which(regresult$coefficients$Model == "formula(fit)")] <- Reduce(paste, deparse(formula(fit)))
    
    # merge data into tracking
    mainregresult$stats <- rbind(mainregresult$stats, regresult$stats)
    mainregresult$changestats <- rbind(mainregresult$changestats, regresult$changestats)
    mainregresult$coefficients <- rbind(mainregresult$coefficients, regresult$coefficients)
  }
  
  # collapse across observations and summarize
  res <- list()
  res$mainregresult <- mainregresult
  res$summary <- Bootstrapcollapseregressions(mainregresult)
  Bootstrapregressionsummary(res$summary)
  
  return(res)
}

Bootstrapcollapseregressions <- function(mainregresult) {
  # collapse across observations
  res <- list()
  tempdbs <- mainregresult$stats
  fun <- function(x){c(m=mean(x, na.rm=TRUE))}
  workingnames <- colnames(tempdbs); workingnames <- workingnames[which(workingnames != 'Model')]
  res$stats <- doBy::summaryBy(tempdbs, formula=list(workingnames, c("Model")), FUN = fun, keep.names = TRUE)
  varsofinterest <- c('r.squared', 'r.squaredadj')
  for (uV in 1:length(varsofinterest)) {
    res$stats[,sprintf('%smin', varsofinterest[uV])] <- NA
    res$stats[,sprintf('%smax', varsofinterest[uV])] <- NA
    uniquemodels <- unique(tempdbs$Model)
    for (uM in 1:length(uniquemodels)) {
      subtempdbs <- tempdbs[which(tempdbs$Model == uniquemodels[uM]),]
      #0.975
      margin <- stats::qt(0.975, df=res$stats$DFd[2])*sd(subtempdbs[,sprintf('%s', varsofinterest[uV])],
                                                                   na.rm=TRUE)/sqrt(res$stats$DFd[2])
      if (is.na(margin)) {
        margin <- 0
      }
      res$stats[which(res$stats$Model == uniquemodels[uM]),
                sprintf('%smin', varsofinterest[uV])] <- mean(subtempdbs[,sprintf('%s', varsofinterest[uV])], na.rm=TRUE) - margin
      res$stats[which(res$stats$Model == uniquemodels[uM]),
                sprintf('%smax', varsofinterest[uV])] <- mean(subtempdbs[,sprintf('%s', varsofinterest[uV])], na.rm=TRUE) + margin
    }
  }
  
  tempdbs <- mainregresult$changestats
  fun <- function(x){c(m=mean(x, na.rm=TRUE))}
  workingnames <- colnames(tempdbs)
  workingnames <- workingnames[which(workingnames != 'Model')]
  workingnames <- workingnames[which(workingnames != 'textoutput')]
  res$changestats <- doBy::summaryBy(tempdbs, formula=list(workingnames, c("Model")), FUN = fun, keep.names = TRUE)
  varsofinterest <- c('r.squared.change')
  for (uV in 1:length(varsofinterest)) {
    res$changestats[,sprintf('%smin', varsofinterest[uV])] <- NA
    res$changestats[,sprintf('%smax', varsofinterest[uV])] <- NA
    uniquemodels <- unique(tempdbs$Model)
    for (uM in 1:length(uniquemodels)) {
      subtempdbs <- tempdbs[which(tempdbs$Model == uniquemodels[uM]),]
      margin <- stats::qt(0.975, df=res$changestats$DFd[2])*sd(subtempdbs[,sprintf('%s', varsofinterest[uV])],
                                                         na.rm=TRUE)/sqrt(res$changestats$DFd[2])
      if (is.na(margin)) {
        margin <- 0
      }
      res$changestats[which(res$changestats$Model == uniquemodels[uM]),
                      sprintf('%smin', varsofinterest[uV])] <- mean(subtempdbs[,sprintf('%s', varsofinterest[uV])], na.rm=TRUE) - margin
      res$changestats[which(res$changestats$Model == uniquemodels[uM]),
                      sprintf('%smax', varsofinterest[uV])] <- mean(subtempdbs[,sprintf('%s', varsofinterest[uV])], na.rm=TRUE) + margin
    }
  }
  
  tempdbs <- mainregresult$coefficients
  fun <- function(x){c(m=mean(x, na.rm=TRUE))}
  workingnames <- colnames(tempdbs)
  workingnames <- workingnames[which(workingnames != 'Model')]
  workingnames <- workingnames[which(workingnames != 'Variable')]
  res$coefficients <- doBy::summaryBy(tempdbs, formula=list(workingnames, c("Model", "Variable")), FUN = fun, keep.names = TRUE)
  res$coefficients$t <- abs(res$coefficients$t)
  res$coefficients$Beta[which(is.na(res$coefficients$Beta))] <- 0
  
  return(res)
}

Bootstrapregressionsummary <- function(res) {
  outputtext <- 'Regression Model Fits:'
  outputtext <- sprintf('%s %s\n', outputtext, res$stats$Model[2])
  tempstr <- Rmimic::fuzzyP(res$changestats$p.value[2])
  # you can get regression p value from F stat using pf(F, DFn, DFd, lower.tail = FALSE) 
  outputtext <- sprintf('%s\tF(%d, %d) = %.1f', outputtext, floor(res$changestats$DFn[2]), ceiling(res$changestats$DFd[2]), res$changestats$F.change[2])
  
  tempstrmod <- ''
  if (as.numeric(tempstr$interpret) <= 0.05) {tempstrmod <- sprintf('%s*', tempstrmod)}
  if (as.numeric(tempstr$interpret) <= 0.01) {tempstrmod <- sprintf('%s*', tempstrmod)}  
  if (as.numeric(tempstr$interpret) <= 0.001) {tempstrmod <- sprintf('%s*', tempstrmod)}  
  outputtext <- sprintf('%s, p %s %s%s', outputtext, tempstr$modifier, tempstr$report, tempstrmod)
  
  outputtext <- sprintf('%s, f\u00b2 = %.2f', outputtext, round(res$changestats$fsquared[2],2))
  outputtext <- sprintf('%s [95%% CI: %.2f to %.2f]', outputtext, round(res$changestats$fsquared.ci.lower[2],2), round(res$changestats$fsquared.ci.upper[2],2))
  outputtext <- sprintf('%s, R\u00b2change = %.2f', outputtext, round(res$changestats$r.squared.change[2],2))
  outputtext <- sprintf('%s [95%% CI: %.2f to %.2f]', outputtext, round(res$changestats$r.squared.changemin[2],2),
                        round(res$changestats$r.squared.changemax[2],2))
  outputtext <- sprintf('%s, R\u00b2adj = %.2f', outputtext, round(res$stats$r.squaredadj[2],2))
  outputtext <- sprintf('%s [95%% CI: %.2f to %.2f]', outputtext, round(res$stats$r.squaredadjmin[2],2),
                        round(res$stats$r.squaredadjmax[2],2))
  outputtext <- sprintf('%s, VIF = %.1f', outputtext, round(res$stats$VIF[2],1))
  cat(outputtext)
  
  outputtext <- '\n\nCoefficients for Model:\n'
  coefdata <- res$coefficients[which(res$coefficients$Model == res$stats$Model[2]),]
  for (cR in 1:nrow(coefdata)) {
    tempstr <- Rmimic::fuzzyP(as.double(coefdata$p.value[cR]))
    coefdata$p.value[cR] <- sprintf('%s %s', tempstr$modifier, tempstr$report)
    
    outputtext <- sprintf('%s\t%s (', outputtext, coefdata$Variable[cR])
    outputtext <- sprintf('%sB = ', outputtext)
    outputtext <- sprintf('%s%.2f', outputtext, round(coefdata$B[cR], digits=2))
    outputtext <- sprintf('%s [95%% CI: %.2f to %.2f]', outputtext, round(coefdata$B.lower.conf.int[cR], digits=2), round(coefdata$B.upper.conf.int[cR], digits=2))
    outputtext <- sprintf('%s, SE B = ', outputtext)
    outputtext <- sprintf('%s%.2f', outputtext, round(coefdata$SE[cR], digits=2))
    if (cR > 1) {
      outputtext <- sprintf('%s, Beta = ', outputtext)
      outputtext <- sprintf('%s%.2f', outputtext, round(coefdata$Beta[cR], digits=2))
    }
    outputtext <- sprintf('%s, t = ', outputtext)
    outputtext <- sprintf('%s%.1f', outputtext, round(coefdata$t[cR], digits=1))
    
    tempstrmod <- ''
    if (as.numeric(tempstr$interpret) <= 0.05) {tempstrmod <- sprintf('%s*', tempstrmod)}
    if (as.numeric(tempstr$interpret) <= 0.01) {tempstrmod <- sprintf('%s*', tempstrmod)}  
    if (as.numeric(tempstr$interpret) <= 0.001) {tempstrmod <- sprintf('%s*', tempstrmod)} 
    outputtext <- sprintf('%s, p %s%s', outputtext, coefdata$p.value[cR], tempstrmod)
    outputtext <- sprintf('%s)\n', outputtext)
  }
  cat(outputtext)
}



