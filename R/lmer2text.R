#' lmer2text
#'
#' @description Output lmerTest results in APA style format with effect sizes and confidence intervals.
#'
#' @param fit lmerTest model fit
#' @param model Parameter to select model type. Currently only ANOVA is available.
#' @param df Parameter to indicate what degrees of freedom approximation should be used. Default is Kenward-Roger. Other option is Shattertwaite.
#' @param numparticipants Number of unique participants in model.
#' @param numfactors Number of fixed and random factors in model.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' 
#' @return A list with:
#' \item{fit}{lmerTest object}
#' \item{ANOVA}{ANOVA table for fixed effects}
#' \item{RandomEffectsANOVA}{ANOVA table for random effects}
#' \item{Rsquared}{model R squared values}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, December 19, 2019
#' 
#' @examples
#' 
#' #Note, use lmerTest::lmer instead of lme4::lmer for inclusion of p values
#' result <- lmer2text(fit, df="Kenward-Roger", numparticipants=30, numfactors=4)
#' 
#' @importFrom stats residuals
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom psychometric CI.Rsq
#' @importFrom lmerTest ranova
#'
#' @export

lmer2text <- function(fit, model=NULL, df=NULL, numparticipants=NULL, numfactors=NULL, confidenceinterval=0.95) {

  spancharacter <- "-"
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
  }  

  result <- list()
  result$fit <- fit
  result$ANOVA <- ''
  result$RandomEffectsANOVA <- ''
  result$Rsquared <- ''
  
  if (is.null(model)) {
    model <- "ANOVA"
  } else {
    model <- "ANOVA"
  }
  
  if (is.null(df)) {
    df = "Kenward-Roger"
  } else {
    if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    } else {
      df = "Kenward-Roger"
    }
  }
  
  # Compute ANOVA from model
  if (toupper(model) == toupper("ANOVA")) {
    as <- NULL
    if (toupper(df) == toupper("Shattertwaite")) {
      as <- data.frame(anova(fit, type = 3))
    } else {
      as <- data.frame(anova(fit, type = 3, ddf = "Kenward-Roger"))
    }
    
    dataframeout <- data.frame(matrix(NA,nrow=nrow(as),ncol=12))
    colnames(dataframeout) <- c("Effect", "DFn", "DFd", "SSn", "SSd", "SSe", "F.value", "p.value", "partialetasquared", "fsquared","fsquared.ci.lower", "fsquared.ci.upper")
    
    dataframeout['Effect'] <- rownames(as)
    dataframeout['DFn'] <- as$NumDF
    dataframeout['DFd'] <- as$DenDF
    dataframeout['SSn'] <- as$Sum.Sq
    dataframeout['SSd'] <- as$Mean.Sq
    dataframeout['SSe'] <- sum(dataframeout$SSd, na.rm=TRUE) + mean(stats::residuals(fit)^2)
    dataframeout['F.value'] <- as$F.value
    dataframeout['p.value'] <- as$Pr..F.
    dataframeout$partialetasquared <- dataframeout$SSn / (dataframeout$SSn + dataframeout$SSe)

    #Calculate f^2 = partialeta^2 / ( 1 - partialeta^2 )
    dataframeout$fsquared <- dataframeout$partialetasquared/(1-dataframeout$partialetasquared)
    for (cT in 1:nrow(dataframeout)) {
      if (!is.na(dataframeout$partialetasquared[cT])) {
        temp <- psychometric::CI.Rsq(dataframeout$partialetasquared[cT], numparticipants, length(unlist(strsplit(as.character(dataframeout$Effect[cT]),"[:]"))), level = confidenceinterval)
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
        rm(temp)
      }
    }
    rm(cT)
    
    result$ANOVA <- dataframeout
    rm(dataframeout)
    
    # place data into text output
    result$ANOVA$textoutput <- NA
    for (cR in 1:nrow(result$ANOVA)) {
      # DF
      pullvalue <- sprintf('%.1f', round(as.numeric(result$ANOVA$DFn[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      temptext <- sprintf('F(%s,', pullvalue)
      
      pullvalue <- sprintf('%.1f', round(as.numeric(result$ANOVA$DFd[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      temptext <- sprintf('%s %s)', temptext, pullvalue)
      # F stat
      pullvalue <- sprintf(' = %.1f', round(as.numeric(result$ANOVA$F.value[cR]), digits = 1))
      if (pullvalue == " = 0.0") {
        pullvalue = " < 0.1"
      }
      temptext <- sprintf('%s%s, p', temptext, pullvalue)
      # P val
      outPvalue <- Rmimic::fuzzyP(as.numeric(result$ANOVA$p.value[cR]))
      temptext <- sprintf('%s %s %s,', temptext, outPvalue$modifier, outPvalue$report)
      rm(outPvalue)
      # effect size
      tempfsq <- sprintf("= %.2f", round(as.numeric(result$ANOVA$fsquared[cR]), digits = 2))
      if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
        tempfsq <- "< 0.01"
      }
      tempfsqlw <- sprintf("%.2f", round(as.numeric(result$ANOVA$fsquared.ci.lower[cR]), digits = 2))
      if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
        tempfsqlw <- "0.0"
      }
      tempfsqup <- sprintf("%.2f", round(as.numeric(result$ANOVA$fsquared.ci.upper[cR]), digits = 2))
      if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
        tempfsqup <- "0.0"
      }
      if (operatingsystem == "Windows") {
        temptext <- sprintf('%s f\u00b2 %s [%2.0f%% CI: %s to %s].', temptext,
                            tempfsq, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      } else {
        temptext <- sprintf('%s f^2 %s [%2.0f%% CI: %s to %s].', temptext,
                            tempfsq, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
      result$ANOVA$textoutput[cR] <- temptext
      rm(tempfsq, tempfsqlw, tempfsqup, temptext)
    }
    rm(cR)
    
    # Compute R^2
    #Nakagawa, S., & Schielzeth, H. (2012). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4, 133-142. doi: 10.1111/j.2041-210x.2012.00261.x
    #Nakagawa, S., Johnson, P. C., & Schielzeth, H. (2017). The coefficient of determination R^2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface, 14(134), 20170213. doi:10.1098/rsif.2017.0213
    dataframeout <- data.frame(matrix(NA,nrow=1,ncol=9))
    colnames(dataframeout) <- c("FixedEffects", "FixedEffects.ci.lower", "FixedEffects.ci.upper", "RandomEffects", "RandomEffects.ci.lower", "RandomEffects.ci.upper", "Model", "Model.ci.lower", "Model.ci.upper")
    
    dataframeout$FixedEffects <- suppressWarnings(MuMIn::r.squaredGLMM(fit)[1])
    temp <- psychometric::CI.Rsq(dataframeout$FixedEffects[1], numparticipants, numfactors, level=confidenceinterval)
    if (temp$LCL[1]<0) {
      temp$LCL[1] <- 0
    } 
    if (temp$UCL[1] > 1) {
      temp$UCL[1] <- Inf
    }
    dataframeout$FixedEffects.ci.lower <- temp$LCL[1]
    dataframeout$FixedEffects.ci.upper <- temp$UCL[1]
    
    dataframeout$Model <- suppressWarnings(MuMIn::r.squaredGLMM(fit)[2])
    temp <- psychometric::CI.Rsq(dataframeout$Model[1], numparticipants, numfactors, level=confidenceinterval)
    if (temp$LCL[1]<0) {
      temp$LCL[1] <- 0
    } 
    if (temp$UCL[1] > 1) {
      temp$UCL[1] <- Inf
    }
    dataframeout$Model.ci.lower <- temp$LCL[1]
    dataframeout$Model.ci.upper <- temp$UCL[1]
    
    dataframeout$RandomEffects <- dataframeout$Model - dataframeout$FixedEffects
    temp <- psychometric::CI.Rsq(dataframeout$RandomEffects[1], numparticipants, numfactors, level=confidenceinterval)
    if (temp$LCL[1]<0) {
      temp$LCL[1] <- 0
    } 
    if (temp$UCL[1] > 1) {
      temp$UCL[1] <- Inf
    }
    dataframeout$RandomEffects.ci.lower <- temp$LCL[1]
    dataframeout$RandomEffects.ci.upper <- temp$UCL[1]
    result$Rsquared <- dataframeout
    
    #https://www.rensvandeschoot.com/tutorials/lme4/
    as <- data.frame(lmerTest::ranova(fit))
    as <- as[2:nrow(as),]
    dataframeout <- data.frame(matrix(NA,nrow=nrow(as),ncol=5))
    colnames(dataframeout) <- c("Effect", "DF", "LogLikelihood", "LRT", "p.value")
    
    dataframeout['Effect'] <- rownames(as)
    dataframeout['DF'] <- as$Df
    dataframeout['LogLikelihood'] <- as$logLik
    dataframeout['LRT'] <- as$LRT
    dataframeout['p.value'] <- as$Pr..Chisq.
    result$RandomEffectsANOVA <- dataframeout
    rm(dataframeout)
    
    # place data into text output
    result$RandomEffectsANOVA$textoutput <- NA
    for (cR in 1:nrow(result$RandomEffectsANOVA)) {
      # DF
      pullvalue <- sprintf('%.1f', round(as.numeric(result$RandomEffectsANOVA$DF[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      temptext <- sprintf('Likelihood Ratio (%s)', pullvalue)
      
      # LRT stat
      pullvalue <- sprintf(' = %.1f', round(as.numeric(result$RandomEffectsANOVA$LRT[cR]), digits = 1))
      if (pullvalue == " = 0.0") {
        pullvalue = " < 0.1"
      }
      temptext <- sprintf('%s%s,', temptext, pullvalue)
      
      # logLik stat
      pullvalue <- sprintf(' = %.1f', round(as.numeric(result$RandomEffectsANOVA$LogLikelihood[cR]), digits = 1))
      if (pullvalue == " = 0.0") {
        pullvalue = " < 0.1"
      }
      temptext <- sprintf('%s log-likelihood%s, p', temptext, pullvalue)
      
      # P val
      outPvalue <- Rmimic::fuzzyP(as.numeric(result$RandomEffectsANOVA$p.value[cR]))
      temptext <- sprintf('%s %s %s,', temptext, outPvalue$modifier, outPvalue$report)
      rm(outPvalue)
      
      result$RandomEffectsANOVA$textoutput[cR] <- temptext
      rm(temptext)
    }
    rm(cR)
  }
  
  return(result)
}