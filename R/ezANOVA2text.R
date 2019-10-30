#' ezANOVA2text
#'
#' @description Output ezANOVA results in APA style format with effect sizes and confidence intervals.
#'
#' @param result ezANOVA output
#' @param numparticipants Number of unique participants in ANOVA
#' @param feffect Parameter to select effect size estimate. Default is Partial Eta Squared. Other option is Generalized Eta Squared.
#' @param sphericity Parameter to select Sphericity correction method. Default is Greenhouse-Geisser. Other option is Huynh-Feldt.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' 
#' @return A list with the same elements as the ezANOVA output:
#' \item{aov}{ANOVA object}
#' \item{ANOVA}{ANOVA table}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 29, 2019
#' 
#' @importFrom psychometric CI.Rsq
#'
#' @export

ezANOVA2text <- function(result, numparticipants=NULL, feffect=NULL, sphericity=NULL, confidenceinterval=0.95, studywiseAlpha=0.05) {
  
  feffect <- NULL
  sphericity <- NULL
  confidenceinterval<-0.95
  studywiseAlpha<-0.05
  numparticipants <- 60
  
  if (is.null(sphericity)) {
    sphericity <- "Greenhouse-Geisser"
  } else {
    if (toupper(sphericity) == toupper("Huynh-Feldt")) {
      sphericity <- "Huynh-Feldt"
    } else {
      sphericity <- "Greenhouse-Geisser"
    }
  }
  if (is.null(feffect)) {
    feffect <- "Partial Eta Squared"
  } else {
    if (toupper(feffect) == toupper("Partial Eta Squared")) {
      feffect <- "Partial Eta Squared"
    } else if (toupper(feffect) == toupper("Generalized Eta Squared")) {
      feffect <- "Generalized Eta Squared"
    } else {
      feffect <- "Partial Eta Squared"
    }
  }
  
  
  spancharacter <- "-"
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
  }  
  
  # Trim off intercept from ANOVA and significance indicator
  result$ANOVA <- result$ANOVA[2:nrow(result$ANOVA),]
  result$ANOVA <- result$ANOVA[,which(colnames(result$ANOVA) != "p<.05")]
  
  # Compute partial eta squared
  result$ANOVA$partialetasquared <- result$ANOVA$SSn/(result$ANOVA$SSn+result$ANOVA$SSd)
  
  # Calculate f^2 = partialeta^2 / ( 1 - partialeta^2 )
  if (toupper(feffect) == toupper("Partial Eta Squared")) {
    result$ANOVA$fsquared <- result$ANOVA$partialetasquared/(1-result$ANOVA$partialetasquared)
  } else {
    result$ANOVA$fsquared <- result$ANOVA$ges/(1-result$ANOVA$ges)
  }
  
  # Compute effect size confidence intervals
  result$ANOVA$fsquared.ci.lower <- NA
  result$ANOVA$fsquared.ci.upper <- NA
  for (cT in 1:nrow(result$ANOVA)) {
    if (toupper(feffect) == toupper("Partial Eta Squared")) {
      tempR2 <- result$ANOVA$partialetasquared[cT]
    } else {
      tempR2 <- result$ANOVA$ges[cT]
    }
    if (!is.na(tempR2)) {
      temp <- psychometric::CI.Rsq(tempR2, numparticipants, length(unlist(strsplit(as.character(result$ANOVA$Effect[cT]),"[:]"))), level = confidenceinterval)
      if (temp$LCL[1]<0) {
        result$ANOVA$fsquared.ci.lower[cT] <- 0
      } else {
        result$ANOVA$fsquared.ci.lower[cT] <- (temp$LCL[1]/(1-temp$LCL[1]))
      }
      if (temp$UCL[1] > 1) {
        result$ANOVA$fsquared.ci.upper[cT] <- Inf
      } else {
        result$ANOVA$fsquared.ci.upper[cT] <- (temp$UCL[1]/(1-temp$UCL[1]))
      }
      rm(temp)
    }
    rm(tempR2)
  }
  rm(cT)
  
  # Check to see if sphericity is violated
  templist <- names(result)
  sphericitytest <- NULL
  sphericitycorrections <- NULL
  if (length(which(templist == "Mauchly's Test for Sphericity"))>0) {
    sphericitytest <- result$`Mauchly's Test for Sphericity`
  }
  if (length(which(templist == "Sphericity Corrections"))>0) {
    sphericitycorrections <- result$`Sphericity Corrections`
  }
  if (!is.null(sphericitytest)) {
    for (cR in 1:nrow(sphericitytest)) {
      if (as.numeric(sphericitytest[cR,"p"]) <= studywiseAlpha) {
        
        cRT <- which(result$ANOVA$Effect == sphericitycorrections$Effect[cR])
        if (toupper(sphericity) == toupper("Greenhouse-Geisser")) {
          result$ANOVA[cRT,"p"] <- sphericitycorrections[cR,"p[GG]"]
          epsilon <- sphericitycorrections[cR,"GGe"]
        } else {
          result$ANOVA[cRT,"p"] <- sphericitycorrections[cR,"p[HF]"]
          epsilon <- sphericitycorrections[cR,"HFe"]
        }
        # update degrees of freedom
        result$ANOVA[cRT,"DFn"] <- sprintf("%.1f", round(as.numeric(result$ANOVA[cRT,"DFn"])/as.numeric(epsilon), digits=1))
        result$ANOVA[cRT,"DFd"] <- sprintf("%.1f", round(as.numeric(result$ANOVA[cRT,"DFd"])/as.numeric(epsilon), digits=1))
        rm(epsilon)
      }
    }
    rm(cR)
  }
  rm(templist, sphericitytest, sphericitycorrections)
  
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
    pullvalue <- sprintf(' = %.1f', round(as.numeric(result$ANOVA$F[cR]), digits = 1))
    if (pullvalue == " = 0.0") {
      pullvalue = " < 0.1"
    }
    temptext <- sprintf('%s%s, p', temptext, pullvalue)
    # P val
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$ANOVA$p[cR]))
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
  
  # populate outputs
  rownames(result$ANOVA) <- 1:nrow(result$ANOVA)
  
  return(result)
}  