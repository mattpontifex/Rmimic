#' RmimicLMcontrast
#'
#' @description Compute SPSS style results for regression analysis with effect size and confidence intervals. This function takes stats::lm fits for a base model and the model of interest and calculates statistics for the model of interest relative to the base model. For logistic regression, pseudo r2 values are reported using Tjurs 2009 approach.
#'
#' @param fit fit from a basic linear model
#' @param altfit fit from a linear model of interest
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#'
#' @return
#' \item{stats}{ANOVA summary table.}
#' \item{changestats}{ANOVA summary table of change statistics.}
#' \item{coefficients}{A summary table of regression coefficients.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 27, 2020
#'
#' @importFrom stats confint extractAIC pf anova coef predict residuals update
#' @importFrom lm.beta lm.beta
#' @importFrom fmsb VIF
#' @importFrom psychometric CI.Rsq
#' @importFrom performance r2_tjur
#' @importFrom utils packageDate
#' @importFrom olsrr ols_test_breusch_pagan
#' 
#'
#' @examples
#'
#'     # Compute Model from mtcars dataset
#'     basefit <- lm(mpg ~ am + wt, data = mtcars)
#'     fit <- lm(mpg ~ am + wt + qsec, data = mtcars)
#'     regresult <- RmimicLMcontrast(basefit, fit, 
#'        confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE)
#'
#' @export

RmimicLMcontrast <- function(fit, altfit, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE) {  
  
  # debug variables
  #fit <- lm(mpg ~ 1, data = mtcars)
  #fit <- lm(mpg ~ am + wt, data = mtcars)
  #altfit <- lm(mpg ~ am + wt, data = mtcars)
  #altfit <- lm(mpg ~ am + wt + qsec, data = mtcars)
  #confidenceinterval<-0.95
  #studywiseAlpha<-0.05
  #verbose<-TRUE
  
  #fit <- glm(vs ~ 1, family = "binomial", data = mtcars)
  #altfit <- glm(vs ~ disp, family = "binomial", data = mtcars)
  
  # check model
  boolsame <- FALSE
  if ((sum(!as.character(fit$call) == as.character(altfit$call))) == 0) {
    boolsame <- TRUE
  }
  
  logisticfit <- FALSE
  if (!is.null(fit$family)) {
    if (toupper(fit$family[1]) == "BINOMIAL") {
      logisticfit <- TRUE
    }
  }
  
  logistic <- FALSE
  if (!is.null(altfit$family)) {
    if (toupper(altfit$family[1]) == "BINOMIAL") {
      logistic <- TRUE
    }
  }
  
  spansize <- 95
  spancharacter <- "-"
  bigspancharacter <- " - "
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
    bigspancharacter <- " _ "
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
    bigspancharacter <- " - "
  }

  # create output structures
  res <- list()
  res$stats <- data.frame(matrix(data=NA,nrow=2,ncol=9))
  names(res$stats) <- c('Model','DFn','DFd','F.value','p.value','r.squared','r.squaredadj','AIC','VIF')
  
  res$modelcheck <- data.frame(matrix(data=NA,nrow=4,ncol=6))
  names(res$modelcheck) <- c('Model', 'Test', 'DF', 'Chisquared', 'p.value', 'decision')
  
  res$changestats <- data.frame(matrix(data=NA,nrow=2,ncol=9))
  names(res$changestats) <- c('Model','r.squared.change','DFn','DFd','F.change','p.value','fsquared', 'fsquared.ci.lower', 'fsquared.ci.upper')
  res$coefficients <- data.frame(matrix(data=NA,nrow=2,ncol=9))
  names(res$coefficients) <- c('Model','Variable','B','SE','Beta','t','p.value','B.lower.conf.int','B.upper.conf.int')
  Coeffcurrentline <- 0

  # populate first model
  mci <- stats::confint(fit, level=confidenceinterval)
  ms <- summary(fit)
  mcoef <- ms$coefficients
  if (logisticfit) {
    mbeta <- exp(stats::coef(fit))
    #Tjur, T. (2009). Coefficients of determination in logistic regression models - A new proposal: The coefficient of discrimination. The American Statistician, 63(4), 366-372.
    ms$r.squared <- performance::r2_tjur(fit)[[1]]
    ms$adj.r.squared <- ms$r.squared
    ms$residuals <- ms$deviance.resid
  } else {
    mbeta <- lm.beta::lm.beta(fit)
  }
  
  ms$f.squared <- ms$r.squared / (1-ms$r.squared)
  temp <- unlist(strsplit(as.character(ms$call), ","))
  res$stats$Model[1] <- temp[2]
  
  if (is.null(ms$fstatistic[3])) {
    if (logisticfit) {
      tempmodel <- anova(stats::update(fit, ~1), fit, test="Chisq")# update here produces null model for comparison
      res$stats$DFn[1] <- abs(tempmodel$Df[2])
      res$stats$DFd[1] <- abs(tempmodel$`Resid. Df`[2])
      res$stats$F.value[1] <- abs(tempmodel$Deviance[2])
      res$stats$p.value[1] <- abs(tempmodel$`Pr(>Chi)`[2])
      if (is.na(abs(tempmodel$`Pr(>Chi)`[2]))) {
        res$stats$p.value[1] <- 1
      }
    } else {
      res$stats$DFn[1] <- 0
      res$stats$DFd[1] <- 0
      res$stats$F.value[1] <- 0
      res$stats$p.value[1] <- 1
    }
  } else {
    res$stats$DFn[1] <- ms$fstatistic[2]
    res$stats$DFd[1] <- ms$fstatistic[3]
    res$stats$F.value[1] <- ms$fstatistic[1]
    res$stats$p.value[1] <- stats::pf(q=ms$fstatistic[1], df1=ms$fstatistic[2], df2=ms$fstatistic[3], lower.tail=FALSE)
  }
  if (is.null(ms$r.squared[1])) {
    res$stats$r.squared[1] <- 0
    res$stats$r.squaredadj[1] <- 0
  } else {
    res$stats$r.squared[1] <- ms$r.squared
    res$stats$r.squaredadj[1] <- ms$adj.r.squared
  }
  res$stats$AIC[1] <- stats::extractAIC(fit)[2]
  res$stats$VIF[1] <- fmsb::VIF(fit)[1]
  
  res$modelcheck$Model[1] <- res$stats$Model[1]
  res$modelcheck$Model[2] <- res$stats$Model[1]
  if (!logisticfit) {
    res$modelcheck$Test[1] <- 'Fitted Values'
    
    tres <- olsrr::ols_test_breusch_pagan(fit) # USES THE FITTED VALUES OF THE MODEL
    res$modelcheck$DF[1] <- 1
    res$modelcheck$Chisquared[1] <- tres$bp
    res$modelcheck$p.value[1] <- tres$p
    #tres <- car::ncvTest(fit)
    #res$modelcheck$DF[1] <- 1
    #res$modelcheck$Chisquared[1] <- tres$ChiSquare
    #res$modelcheck$p.value[1] <- tres$p
    
    res$modelcheck$Test[2] <- 'Indep Variables'
    
    tres <- olsrr::ols_test_breusch_pagan(fit, rhs = TRUE) # USES THE INDEPENDENT VARIABLES OF THE MODEL
    res$modelcheck$DF[2] <- length(tres$preds)
    res$modelcheck$Chisquared[2] <- tres$bp
    res$modelcheck$p.value[2] <- tres$p
    #tres <- lmtest::bptest(fit, studentize=TRUE)
    #res$modelcheck$DF[2] <- tres$parameter[[1]]
    #res$modelcheck$Chisquared[2] <- tres$statistic[[1]]
    #res$modelcheck$p.value[2] <- tres$p.value[[1]]
    
  }
  
  res$changestats$Model[1] <- res$stats$Model[1]
  res$changestats$r.squared.change[1] <- res$stats$r.squared[1]
  res$changestats$DFn[1] <- res$stats$DFn[1]
  res$changestats$DFd[1] <- res$stats$DFd[1]
  res$changestats$F.change[1] <- res$stats$F.value[1]
  res$changestats$p.value[1] <- res$stats$p.value[1]
  
  if (is.null(ms$r.squared[1])) {
    res$changestats$fsquared[1] <- 0
  } else {
    res$changestats$fsquared[1] <- ms$f.squared
  }
  if (!is.null(ms$r.squared[1])) {
    temporig <- psychometric::CI.Rsq(ms$r.squared, length(ms$residuals), length(trimws(unlist(strsplit(as.character(ms$terms[[3]]),"[+]"))))-1, level = confidenceinterval)
    if ((temporig$LCL[1]<0) | (is.na(temporig$LCL[1]))) {
      res$changestats$fsquared.ci.lower[1] <- 0
    } else {
      res$changestats$fsquared.ci.lower[1] <- (temporig$LCL[1]/(1-temporig$LCL[1]))
    }
    if ((temporig$UCL[1] > 1) | (is.na(temporig$UCL[1]))) {
      res$changestats$fsquared.ci.upper[1] <- Inf
    } else {
      res$changestats$fsquared.ci.upper[1] <- (temporig$UCL[1]/(1-temporig$UCL[1]))
    }
  }
  
  for (cV in 1:nrow(mcoef)) {
    res$coefficients[Coeffcurrentline + cV, 1] <- res$stats$Model[1]
    res$coefficients$Variable[Coeffcurrentline + cV] <- row.names(mcoef)[cV]
    res$coefficients$B[Coeffcurrentline + cV] <- mcoef[cV,1]
    res$coefficients$SE[Coeffcurrentline + cV] <- mcoef[cV,2]
    if (logisticfit) {
      res$coefficients$Beta[Coeffcurrentline + cV] <- mbeta[cV]
    } else {
      res$coefficients$Beta[Coeffcurrentline + cV] <- mbeta$standardized.coefficients[cV]
    }
    res$coefficients$t[Coeffcurrentline + cV] <- mcoef[cV,3]
    res$coefficients$p.value[Coeffcurrentline + cV] <- mcoef[cV,4]
    if (logisticfit) {
      if (nrow(mcoef) == 1) {
        res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[1]
        res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[2]
      } else {
        res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[cV,1]
        res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[cV,2]
      }
    } else {
      res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[cV,1]
      res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[cV,2]
    }
  }
  Coeffcurrentline <- Coeffcurrentline + nrow(mcoef)
  
  # populate second model
  mci <- stats::confint(altfit, level=confidenceinterval)
  msalt <- summary(altfit)
  mcoef <- msalt$coefficients
  if (logistic) {
    mbeta <- exp(stats::coef(altfit))
    #Tjur, T. (2009). Coefficients of determination in logistic regression models - A new proposal: The coefficient of discrimination. The American Statistician, 63(4), 366-372.
    msalt$r.squared <- performance::r2_tjur(altfit)[[1]]
    msalt$adj.r.squared <- msalt$r.squared
    msalt$residuals <- msalt$deviance.resid
  } else {
    mbeta <- lm.beta::lm.beta(altfit)
  }
  msalt$f.squared <- msalt$r.squared / (1-msalt$r.squared)
  temp <- unlist(strsplit(as.character(msalt$call), ","))
  res$stats$Model[2] <- temp[2]
  if (is.null(msalt$fstatistic[3])) {
    if (logistic) {
      tempmodel <- anova(update(altfit, ~1), altfit, test="Chisq")# update here produces null model for comparison
      res$stats$DFn[2] <- abs(tempmodel$Df[2])
      res$stats$DFd[2] <- abs(tempmodel$`Resid. Df`[2])
      res$stats$F.value[2] <- abs(tempmodel$Deviance[2])
      res$stats$p.value[2] <- abs(tempmodel$`Pr(>Chi)`[2])
      if (is.na(abs(tempmodel$`Pr(>Chi)`[2]))) {
        res$stats$p.value[2] <- 1
      }
    } else {
      res$stats$DFn[2] <- 0
      res$stats$DFd[2] <- 0
      res$stats$F.value[2] <- 0
      res$stats$p.value[2] <- 1
    }
  } else {
    res$stats$DFn[2] <- msalt$fstatistic[2]
    res$stats$DFd[2] <- msalt$fstatistic[3]
    res$stats$F.value[2] <- msalt$fstatistic[1]
    res$stats$p.value[2] <- stats::pf(q=msalt$fstatistic[1], df1=msalt$fstatistic[2], df2=msalt$fstatistic[3], lower.tail=FALSE)
  }
  res$stats$r.squared[2] <- msalt$r.squared
  res$stats$r.squaredadj[2] <- msalt$adj.r.squared
  res$stats$AIC[2] <- stats::extractAIC(altfit)[2]
  res$stats$VIF[2] <- fmsb::VIF(altfit)[1]
  
  res$modelcheck$Model[3] <- res$stats$Model[2]
  res$modelcheck$Model[4] <- res$stats$Model[2]
  if (!logisticfit) {
    res$modelcheck$Test[3] <- 'Fitted Values'
    
    tres <- olsrr::ols_test_breusch_pagan(altfit) # USES THE FITTED VALUES OF THE MODEL
    res$modelcheck$DF[3] <- 1
    res$modelcheck$Chisquared[3] <- tres$bp
    res$modelcheck$p.value[3] <- tres$p
    #tres <- car::ncvTest(fit)
    #res$modelcheck$DF[3] <- 1
    #res$modelcheck$Chisquared[3] <- tres$ChiSquare
    #res$modelcheck$p.value[3] <- tres$p
    
    res$modelcheck$Test[4] <- 'Indep Variables'
    
    tres <- olsrr::ols_test_breusch_pagan(altfit, rhs = TRUE) # USES THE INDEPENDENT VARIABLES OF THE MODEL
    res$modelcheck$DF[4] <- length(tres$preds)
    res$modelcheck$Chisquared[4] <- tres$bp
    res$modelcheck$p.value[4] <- tres$p
    #tres <- lmtest::bptest(fit, studentize=TRUE)
    #res$modelcheck$DF[4] <- tres$parameter[[1]]
    #res$modelcheck$Chisquared[4] <- tres$statistic[[1]]
    #res$modelcheck$p.value[4] <- tres$p.value[[1]]
    
  }
  
  
  changeresult <- stats::anova(fit, altfit)
  if (logistic) {
    changeresult <- anova(fit, altfit, test="Chisq")
  }
  res$changestats$Model[2] <- res$stats$Model[2]
  res$changestats$r.squared.change[2] <- res$stats$r.squared[2] - res$stats$r.squared[1]
  res$changestats$DFn[2] <- res$stats$DFn[2]
  res$changestats$DFd[2] <- res$stats$DFd[2]
  if (logistic) {
    res$changestats$F.change[2] <- abs(changeresult$Deviance[2])
    res$changestats$p.value[2] <- abs(changeresult$`Pr(>Chi)`[2])
  } else {
    res$changestats$F.change[2] <- changeresult$F[-1]
    res$changestats$p.value[2] <- changeresult$`Pr(>F)`[-1]
  }
  
  if (is.na(res$changestats$F.change[1])) {
    res$changestats$F.change[1] <- 0.0
  }
  if (is.na(res$changestats$F.change[2])) {
    res$changestats$F.change[2] <- 0.0
  }
  if (is.na(res$changestats$p.value[1])) {
    res$changestats$p.value[1] <- 0.99999
  }
  if (is.na(res$changestats$p.value[2])) {
    res$changestats$p.value[2] <- 0.99999
  }
  
  res$changestats$fsquared[2] <- res$changestats$r.squared.change[2]/(1-res$stats$r.squared[1])
  temp <- psychometric::CI.Rsq(res$changestats$r.squared.change[2], length(msalt$residuals), length(trimws(unlist(strsplit(as.character(msalt$terms[[3]]),"[+]"))))-1, level = confidenceinterval)
  if ((temp$LCL[1]<0) | (is.na(temp$LCL[1]))) {
    res$changestats$fsquared.ci.lower[2] <- 0
  } else {
    res$changestats$fsquared.ci.lower[2] <- (temp$LCL[1]/(1-temporig$LCL[1]))
  }
  if ((temp$UCL[1] > 1) | (is.na(temp$UCL[1]))) {
    res$changestats$fsquared.ci.upper[2] <- Inf
  } else {
    res$changestats$fsquared.ci.upper[2] <- (temp$UCL[1]/(1-temporig$UCL[1]))
  }
  if (boolsame == FALSE) {
    for (cV in 1:nrow(mcoef)) {
      res$coefficients[Coeffcurrentline + cV, 1] <- res$stats$Model[2]
      res$coefficients$Variable[Coeffcurrentline + cV] <- row.names(mcoef)[cV]
      res$coefficients$B[Coeffcurrentline + cV] <- mcoef[cV,1]
      res$coefficients$SE[Coeffcurrentline + cV] <- mcoef[cV,2]
      if (logistic) {
        res$coefficients$Beta[Coeffcurrentline + cV] <- mbeta[cV]
      } else {
        res$coefficients$Beta[Coeffcurrentline + cV] <- mbeta$standardized.coefficients[cV]
      }
      res$coefficients$t[Coeffcurrentline + cV] <- mcoef[cV,3]
      res$coefficients$p.value[Coeffcurrentline + cV] <- mcoef[cV,4]
      if (logistic) {
        if (nrow(mcoef) == 1) {
          res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[1]
          res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[2]
        } else {
          res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[cV,1]
          res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[cV,2]
        }
      } else {
        res$coefficients$B.lower.conf.int[Coeffcurrentline + cV] <- mci[cV,1]
        res$coefficients$B.upper.conf.int[Coeffcurrentline + cV] <- mci[cV,2]
      }
      
    }
    Coeffcurrentline <- Coeffcurrentline + nrow(mcoef)
  }
  
  # interpret modelcheck
  for (cR in 1:nrow(res$modelcheck)) {
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$modelcheck$p.value[cR]))
    if (outPvalue$interpret <= studywiseAlpha) {
      res$modelcheck$decision[cR] <- 'Unacceptable.'
    } else {
      res$modelcheck$decision[cR] <- 'Acceptable.'
    }
  }
  
  # place data into text output
  res$changestats$textoutput <- NA
  cR <- 1
  for (cR in 1:nrow(res$changestats)) {
    # DF
    pullvalue <- sprintf('%.1f', round(as.numeric(res$changestats$DFn[cR]), digits = 1))
    if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
      if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
        pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
      }
    }
    if (logistic) {
      if (cR == 1) {
        temptext <- sprintf('Chi-square(%s,', pullvalue)
      } else {
        temptext <- sprintf('Chi-square change(%s,', pullvalue)
      }
    } else {
      if (cR == 1) {
        temptext <- sprintf('F(%s,', pullvalue)
      } else {
        temptext <- sprintf('Fchange(%s,', pullvalue)
      }
    }
    pullvalue <- sprintf('%.1f', round(as.numeric(res$changestats$DFd[cR]), digits = 1))
    if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
      if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
        pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
      }
    }
    temptext <- sprintf('%s %s)', temptext, pullvalue)
    # F stat
    pullvalue <- sprintf(' = %.1f', round(as.numeric(res$changestats$F.change[cR]), digits = 1))
    if (pullvalue == " = 0.0") {
      pullvalue = " < 0.1"
    }
    temptext <- sprintf('%s%s, p', temptext, pullvalue)
    # P val
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$changestats$p.value[cR]))
    temptext <- sprintf('%s %s %s,', temptext, outPvalue$modifier, outPvalue$report)
    rm(outPvalue)
    # effect size
    tempfsq <- sprintf("= %.2f", round(as.numeric(res$changestats$fsquared[cR]), digits = 2))
    if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
      tempfsq <- "< 0.01"
    }
    tempfsqlw <- sprintf("%.2f", round(as.numeric(res$changestats$fsquared.ci.lower[cR]), digits = 2))
    if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
      tempfsqlw <- "0.0"
    }
    tempfsqup <- sprintf("%.2f", round(as.numeric(res$changestats$fsquared.ci.upper[cR]), digits = 2))
    if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
      tempfsqup <- "0.0"
    }
    if (operatingsystem == "Windows") {
      temptext <- sprintf('%s f\u00b2 %s [%2.0f%% CI: %s to %s]', temptext,
                          tempfsq, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
    } else {
      temptext <- sprintf('%s f^2 %s [%2.0f%% CI: %s to %s]', temptext,
                          tempfsq, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
    }
    res$changestats$textoutput[cR] <- temptext
    rm(tempfsq, tempfsqlw, tempfsqup, temptext)
  }
  
  if (verbose == TRUE) {
    
    temptext <- "Regression Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- "Analysis were conducted using the"
    outstring <- sprintf('%s stats (R Core Team, %s)', outstring, strsplit(as.character(utils::packageDate("stats")),"-")[[1]][1])
    outstring <- sprintf('%s, fmsb (Nakazawa, %s)', outstring, strsplit(as.character(utils::packageDate("fmsb")),"-")[[1]][1])
    outstring <- sprintf('%s, psychometric (Fletcher, %s)', outstring, strsplit(as.character(utils::packageDate("psychometric")),"-")[[1]][1])
    
    if (logistic) {
      outstring <- sprintf('%s, performance (Lüdecke, Makowski, Waggoner, & Patil, %s)', outstring, strsplit(as.character(utils::packageDate("performance")),"-")[[1]][1])
    } else {
      outstring <- sprintf('%s, lm.beta (Behrendt, %s)', outstring, strsplit(as.character(utils::packageDate("lm.beta")),"-")[[1]][1])
    }
    outstring <- sprintf('%s, and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    # show model
    outstring <- "Regression Model Fit"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    outputdataframe <- res$stats
    #c('Model','DFn','DFd','F.value','p.value','r.squared','r.squaredadj','AIC', 'VIF')
    for (cC in 2:3) {
      outputdataframe[,cC] <- sprintf('%0.1f', round(outputdataframe[,cC], digits=1))
    }
    outputdataframe[,4] <- sprintf('%0.1f', round(outputdataframe[,4], digits=1))
    for (cC in 6:7) {
      outputdataframe[,cC] <- sprintf('%0.3f', round(outputdataframe[,cC], digits=3))
    }
    for (cC in 8:9) {
      outputdataframe[,cC] <- sprintf('%0.1f', round(outputdataframe[,cC], digits=1))
    }
    
    # manage specific values
    for (cR in 1:nrow(outputdataframe)) {
      # check F
      if (outputdataframe$F.value[cR] == "0.0") {
        outputdataframe$F.value[cR] <- "< 0.1"
      }
      
      # DF
      pullvalue <- outputdataframe$DFn[cR]
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      outstring = sprintf("%s", pullvalue)
      pullvalue <- outputdataframe$DFd[cR]
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      outstring = sprintf("%s, %s", outstring, pullvalue)
      outputdataframe$DFn[cR] <- outstring
      
      outPvalue <- Rmimic::fuzzyP(as.double(outputdataframe[cR,5]))
      if (outPvalue$modifier == "=") {
        pullvalue <- outPvalue$report
      } else {
        pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      if (outPvalue$interpret <= studywiseAlpha) {
        pullvalue <- sprintf('%s%s', pullvalue, "**")
      }
      outputdataframe[cR,5] <- pullvalue
      rm(outPvalue)
    }
    
    outputdataframe <- outputdataframe[,c('Model','DFn','F.value','p.value','r.squared','r.squaredadj','AIC', 'VIF')]
    if (logistic) {
      names(outputdataframe) <- c('Model','df','Chisq','p','r.squared','r.squaredadj','AIC', 'VIF')
    } else {
      names(outputdataframe) <- c('Model','df','F','p','r.squared','r.squaredadj','AIC', 'VIF')
    }
    sepgap <- data.frame(matrix(floor(spansize/ncol(outputdataframe)), nrow=1, ncol=ncol(outputdataframe)))
    sepgap[1,1] = sepgap[1,1] + 3
    sepgap[1,2] = sepgap[1,2] + 2
    sepgap[1,4] = sepgap[1,4] + 2
    sepgap[1,5] = sepgap[1,5] + 2
    sepgap[1,6] = sepgap[1,6] + 4
    sepgap[1,8] = sepgap[1,8] - 6
    Rmimic::table2console(outputdataframe, sepgap=sepgap, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
    rm(outputdataframe)
    
    if (!logistic) {
      outstring <- "Heteroscedasticity Model Check"
      cat(sprintf("\n"))
      Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
      rm(outstring)
      outputdataframe <- res$modelcheck
      #names(res$modelcheck) <- c('Model', 'Test', 'DF', 'Chisquared', 'p.value', 'decision')
      names(outputdataframe) <- c('Model', 'Using', 'df', 'Chisq', 'p', 'decision')
      outputdataframe[,4] <- sprintf('%0.1f', round(outputdataframe[,4], digits=1))
      for (cR in 1:nrow(outputdataframe)) {
        outPvalue <- Rmimic::fuzzyP(as.double(outputdataframe[cR,5]))
        if (outPvalue$modifier == "=") {
          pullvalue <- outPvalue$report
        } else {
          pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
        }
        if (outPvalue$interpret <= studywiseAlpha) {
          pullvalue <- sprintf('%s%s', pullvalue, "**")
        }
        outputdataframe[cR,5] <- pullvalue
      }
      sepgap <- data.frame(matrix(floor(spansize/ncol(outputdataframe)), nrow=1, ncol=ncol(outputdataframe)))
      sepgap[1,1] = sepgap[1,1] - 3
      sepgap[1,2] = sepgap[1,2] + 9
      sepgap[1,3] = sepgap[1,3] - 4
      sepgap[1,4] = sepgap[1,4] - 4
      sepgap[1,5] = sepgap[1,5] - 3
      sepgap[1,6] = sepgap[1,6] + 9
      Rmimic::table2console(outputdataframe, sepgap=sepgap, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
      rm(outputdataframe)
    }
    
    outstring <- "Change Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    outputdataframe <- res$changestats[,c('Model','r.squared.change','DFn','DFd','F.change','p.value','fsquared', 'fsquared.ci.lower', 'fsquared.ci.upper')]
    outputdataframe[,2] <- sprintf('%0.3f', round(outputdataframe[,2], digits=3))
    for (cC in 3:5) {
      outputdataframe[,cC] <- sprintf('%0.1f', round(outputdataframe[,cC], digits=1))
    }
    for (cC in 7:9) {
      outputdataframe[,cC] <- sprintf('%0.2f', round(outputdataframe[,cC], digits=2))
    }
    # manage specific values
    for (cR in 1:nrow(outputdataframe)) {
      # check F
      if (outputdataframe$F.change[cR] == "0.0") {
        outputdataframe$F.change[cR] <- "< 0.1"
      }
      
      # DF
      pullvalue <- outputdataframe$DFn[cR]
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      outstring = sprintf("%s", pullvalue)
      pullvalue <- outputdataframe$DFd[cR]
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      outstring = sprintf("%s, %s", outstring, pullvalue)
      outputdataframe$DFn[cR] <- outstring
      
      # fsquared
      outstring = sprintf("%s", outputdataframe$fsquared[cR])
      pullvalue <- outputdataframe$fsquared.ci.lower[cR]
      outstring = sprintf("%s [%s", outstring, pullvalue)
      pullvalue <- outputdataframe$fsquared.ci.upper[cR]
      outstring = sprintf("%s to %s]", outstring, pullvalue)
      outputdataframe$fsquared[cR] <- outstring
      
      outPvalue <- Rmimic::fuzzyP(as.double(outputdataframe$p.value[cR]))
      if (outPvalue$modifier == "=") {
        pullvalue <- outPvalue$report
      } else {
        pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      if (outPvalue$interpret <= studywiseAlpha) {
        pullvalue <- sprintf('%s%s', pullvalue, "**")
      }
      outputdataframe$p.value[cR] <- pullvalue
      rm(outPvalue)
    }
    outputdataframe <- outputdataframe[,c('Model','DFn','F.change','p.value','r.squared.change','fsquared')]
    if (logistic) {
      names(outputdataframe) <- c('Model','df','Chisq','p','r.squared','fsquared')
    } else {
      names(outputdataframe) <- c('Model','df','F.change','p','r.squared','fsquared')
    }
    sepgap <- data.frame(matrix(floor(spansize/ncol(outputdataframe)), nrow=1, ncol=ncol(outputdataframe)))
    sepgap[1,1] = sepgap[1,1] - 2
    sepgap[1,2] = sepgap[1,2] - 3
    sepgap[1,3] = sepgap[1,3] - 2
    sepgap[1,4] = sepgap[1,4] - 0
    sepgap[1,5] = sepgap[1,5] - 0
    sepgap[1,6] = sepgap[1,6] + 7
    Rmimic::table2console(outputdataframe, sepgap=sepgap, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
    rm(outputdataframe)
    
    
    outputdataframe <- res$coefficients[,c('Model','Variable','B','B.lower.conf.int','B.upper.conf.int','SE','Beta','t','p.value')]
    for (cC in 3:7) {
      outputdataframe[,cC] <- sprintf('%0.2f', round(outputdataframe[,cC], digits=2))
    }
    outputdataframe[,8] <- sprintf('%0.1f', round(outputdataframe[,8], digits=1))
    # manage specific values
    for (cR in 1:nrow(outputdataframe)) {
      # check t
      if (outputdataframe$t[cR] == "0.0") {
        outputdataframe$t[cR] <- "< 0.1"
      }
      
      # B
      outstring = sprintf("%s", outputdataframe$B[cR])
      pullvalue <- outputdataframe$B.lower.conf.int[cR]
      outstring = sprintf("%s [%s", outstring, pullvalue)
      pullvalue <- outputdataframe$B.upper.conf.int[cR]
      outstring = sprintf("%s to %s]", outstring, pullvalue)
      outputdataframe$B[cR] <- outstring
      
      outPvalue <- Rmimic::fuzzyP(as.double(outputdataframe$p.value[cR]))
      if (outPvalue$modifier == "=") {
        pullvalue <- outPvalue$report
      } else {
        pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      if (outPvalue$interpret <= studywiseAlpha) {
        pullvalue <- sprintf('%s%s', pullvalue, "**")
      }
      outputdataframe$p.value[cR] <- pullvalue
      rm(outPvalue)
    }
    outputdataframe <- outputdataframe[,c('Model','Variable','B','SE','Beta','t','p.value')]
    names(outputdataframe) <- c('Model','Variable','B','SE','Beta','t','p')
    
    modelcalls <- unique(as.character(outputdataframe$Model))
    for (cR in 1:length(modelcalls)) {
      cat(sprintf("\n"))
      outstring <- sprintf("Coefficients for Model: %s", modelcalls[cR])
      Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
      rm(outstring)
      
      suboutputdataframe <- outputdataframe[which(outputdataframe$Model == modelcalls[cR]),c('Variable','B','SE','Beta','t','p')]
      if (logistic) {
        names(suboutputdataframe) <- c('Variable','B','SE','OR','z','p')
      }
      sepgap <- data.frame(matrix(floor(spansize/ncol(suboutputdataframe)), nrow=1, ncol=ncol(suboutputdataframe)))
      Rmimic::table2console(suboutputdataframe, sepgap=NULL, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
      rm(suboutputdataframe)
    }
    
    
    outstring <- "Interpretation with Test Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    outputdataframe <- res$coefficients[,c('Model','Variable','B','B.lower.conf.int','B.upper.conf.int','SE','Beta','t','p.value')]
    
    outtextstring <- sprintf("Step 1: Model: %s", modelcalls[1])
    Rmimic::typewriter(outtextstring, tabs=0, spaces=0, characters=floor(spansize *.90))
    
    suboutputdataframe <- outputdataframe[which(outputdataframe$Model == modelcalls[1]),]
    outtextstring <- sprintf("Regression analysis indicated that")
    
    if (!is.na(suboutputdataframe$Variable[2])) {
      for (cCoef in 2:nrow(suboutputdataframe)) {
        outtextstring <- sprintf("%s %s", outtextstring, suboutputdataframe$Variable[cCoef]) # add variable label
        outtextstring <- sprintf("%s (B = %1.2f", outtextstring, suboutputdataframe$B[cCoef]) # add B
        outtextstring <- sprintf("%s [%2.0f%% CI: %1.2f to %1.2f],", outtextstring, confidenceinterval*100, suboutputdataframe$B.lower.conf.int[cCoef], suboutputdataframe$B.upper.conf.int[cCoef]) # add B
        outtextstring <- sprintf("%s SE B = %1.2f,", outtextstring, suboutputdataframe$SE[cCoef]) # add SE B
        if (logistic) {
          outtextstring <- sprintf("%s Odds Ratio = %1.2f)", outtextstring, suboutputdataframe$Beta[cCoef]) # add OR
        } else {
          outtextstring <- sprintf("%s Beta = %1.2f)", outtextstring, suboutputdataframe$Beta[cCoef]) # add Beta
        }
        if ((nrow(suboutputdataframe)-cCoef) > 0) {
          if ((nrow(suboutputdataframe)-cCoef) == 1) {
            if (nrow(suboutputdataframe) > 3) {
              outtextstring <- sprintf("%s, and", outtextstring)
            } else {
              outtextstring <- sprintf("%s and", outtextstring)
            }
          } else {
            outtextstring <- sprintf("%s, ", outtextstring)
          }
        }
      }
    } else {
      outtextstring <- sprintf("%s the inclusion of a constant", outtextstring)
    }    
      
    temppval <- Rmimic::fuzzyP(as.numeric(res$stats$p.value[1]))
    outval <- paste(temppval$modifier,temppval$interpret,sep = " ")
    if (temppval$interpret <= studywiseAlpha) {
      outtextstring <- sprintf("%s explained a statistically significant", outtextstring)
    } else {
      if (nrow(suboutputdataframe) > 2) {
        outtextstring <- sprintf("%s do not explain a statistically significant", outtextstring)
      } else {
        outtextstring <- sprintf("%s did not explain a statistically significant", outtextstring)
      }
    }
    outtextstring <- sprintf("%s (%s),", outtextstring, res$changestats$textoutput[1])
    outtextstring <- sprintf("%s amount of variance in %s", outtextstring, trimws(strsplit(modelcalls[1], split = '~')[[1]][1]))
    outtextstring <- sprintf("%s (R^2adj = %.2f).", outtextstring, res$stats$r.squaredadj[1])
    
    if (!is.na(res$stats$VIF[1])) {
      if (res$stats$VIF[1] > 2.5) {
        outtextstring <- sprintf("%s [Warning: VIF is %.1f, suggesting a", outtextstring, res$stats$VIF[1])
        if (res$stats$VIF[1] > 2.5) {
          degree <- "slight"
        }
        if (res$stats$VIF[1] > 5) {
          degree <- "moderate"
        }
        if (res$stats$VIF[1] > 10) {
          degree <- "severe"
        }
        outtextstring <- sprintf("%s %s", outtextstring, degree)
        outtextstring <- sprintf("%s multicollinearity issue].", outtextstring)
      }
    }
    outtextstring <- sprintf("%s\n", outtextstring)
    Rmimic::typewriter(outtextstring, tabs=1, spaces=0, characters=floor(spansize *.90))
    rm(outtextstring)
    
    if (boolsame == FALSE) {
      outtextstring <- sprintf("Step 2: Model: %s", modelcalls[2])
      Rmimic::typewriter(outtextstring, tabs=0, spaces=0, characters=floor(spansize *.90))
      outtextstring <- sprintf("Hierarchical regression analysis indicated that")
      suboutputdataframe <- outputdataframe[which(outputdataframe$Model == modelcalls[2]),]
      for (cCoef in 2:nrow(suboutputdataframe)) {
        outtextstring <- sprintf("%s %s", outtextstring, suboutputdataframe$Variable[cCoef]) # add variable label
        outtextstring <- sprintf("%s (B = %1.2f", outtextstring, suboutputdataframe$B[cCoef]) # add B
        outtextstring <- sprintf("%s [%2.0f%% CI: %1.2f to %1.2f],", outtextstring, confidenceinterval*100, suboutputdataframe$B.lower.conf.int[cCoef], suboutputdataframe$B.upper.conf.int[cCoef]) # add B
        outtextstring <- sprintf("%s SE B = %1.2f,", outtextstring, suboutputdataframe$SE[cCoef]) # add SE B
        if (logistic) {
          outtextstring <- sprintf("%s Odds Ratio = %1.2f)", outtextstring, suboutputdataframe$Beta[cCoef]) # add OR
        } else {
          outtextstring <- sprintf("%s Beta = %1.2f)", outtextstring, suboutputdataframe$Beta[cCoef]) # add Beta
        }
        if ((nrow(suboutputdataframe)-cCoef) > 0) {
          if ((nrow(suboutputdataframe)-cCoef) == 1) {
            if (nrow(suboutputdataframe) > 3) {
              outtextstring <- sprintf("%s, and", outtextstring)
            } else {
              outtextstring <- sprintf("%s and", outtextstring)
            }
          } else {
            outtextstring <- sprintf("%s, ", outtextstring)
          }
        }
      }
      
      temppval <- Rmimic::fuzzyP(as.numeric(res$changestats$p.value[2]))
      outval <- paste(temppval$modifier,temppval$interpret,sep = " ")
      if (temppval$interpret <= studywiseAlpha) {
        outtextstring <- sprintf("%s explained a statistically significant", outtextstring)
      } else {
        outtextstring <- sprintf("%s failed to explain a statistically significant", outtextstring)
      }
      outtextstring <- sprintf("%s (%s),", outtextstring, res$changestats$textoutput[2])
      outtextstring <- sprintf("%s change in variance in %s", outtextstring, trimws(strsplit(modelcalls[2], split = '~')[[1]][1]))
      outtextstring <- sprintf("%s (R^2change = %.2f, R^2adj = %.2f).", outtextstring, res$changestats$r.squared.change[2], res$stats$r.squaredadj[2])
      
      if (!is.na(res$stats$VIF[1])) {
        if (res$stats$VIF[2] > 2.5) {
          outtextstring <- sprintf("%s [Warning: VIF is %.1f, suggesting a", outtextstring, res$stats$VIF[2])
          if (res$stats$VIF[2] > 2.5) {
            degree <- "slight"
          }
          if (res$stats$VIF[2] > 5) {
            degree <- "moderate"
          }
          if (res$stats$VIF[2] > 10) {
            degree <- "severe"
          }
          outtextstring <- sprintf("%s %s", outtextstring, degree)
          outtextstring <- sprintf("%s multicollinearity issue].", outtextstring)
        }
      }
      
      outtextstring <- sprintf("%s\n", outtextstring)
      Rmimic::typewriter(outtextstring, tabs=1, spaces=0, characters=floor(spansize *.90))
    }
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  return(res)
}