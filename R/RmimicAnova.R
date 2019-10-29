#' RmimicAnova
#'
#' @description Compute SPSS style univariate ANOVA with effect size and confidence intervals using the ezANOVA function. Main effects and interactions are automatically decomposed using the specified post-hoc corrections. Note: Data should be in long format.
#'
#' @usage
#' result <- RmimicAnova(
#'    data = database,
#'    dependentvariable = "dv",
#'    subjectid = "id",
#'    between = "betweenfactor",
#'    within = "withinfactor")
#'
#' @param data Database containing data
#' @param dependentvariable Dependent Variable label
#' @param subjectid Subject ID label.
#' @param between Between subjects column labels.
#' @param within Within subjects column labels.
#' @param sphericity Parameter to select Sphericity correction method. Default is Greenhouse-Geisser. Other option is Huynh-Feldt.
#' @param posthoc Parameter to indicate what post-hoc comparisons should be performed. Default is Bonferroni. Other options are Holm-Bonferroni, Scheffe, Sidak, Tukey, or False Discovery Rate Control.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param strict Boolean operator for if reported degrees of freedom should be corrected. Default is TRUE.
#' @param printoutput Boolean operator for if the statistics should be printed to the console. Default is TRUE.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#' @param output Parameter to specify if only the ANOVA table should be outputted ("Short") or if descriptive statistics should be outputted as well ("Full").
#' @param planned Parameter to specify an effect to show the post-hoc comparisons even if they are not significant.
#'
#' @return
#' \item{stats}{ANOVA summary table.}
#' \item{aov}{An aov object corresponding to the requested ANOVA.}
#' \item{posthocANOVA}{ANOVA summary table for posthoc decompositions.}
#' \item{posthocttest}{A t test summary table for posthoc decompositions.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 21, 2019
#'
#' @importFrom ez ezANOVA
#' @importFrom psychometric CI.Rsq
#' @importFrom doBy summaryBy
#'
#' @export

RmimicAnova <- function(data, dependentvariable=NULL, subjectid=NULL, between=NULL, within=NULL, nonparametric=FALSE, feffect=NULL, sphericity=NULL, posthoc=NULL, planned=NULL, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE, currentlevel=NULL, verbosedescriptives=TRUE) {

  # revise to incorporate posthoc adjustments
  # revise to incorporate data screening to provide useful error information.
  
  options(contrasts = c("contr.sum", "contr.poly"))
  oldw <- getOption("warn")
  options(warn = -1) # ezANOVA likes to warn about lots of crap

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
  
  if (!is.null(posthoc)) {
    if (toupper(posthoc) == toupper("Bonferroni")) {
      posthoc <- "Bonferroni"
    } else if (toupper(posthoc) == toupper("Sidak")) {
      posthoc <- "Sidak"
    } else if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
      posthoc <- "Holm-Bonferroni"
    } else if (posthoc == FALSE) {
      posthoc <- NULL
    } else {
      posthoc <- "Bonferroni"
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
  
  # Assess what got fed into the function
  betweenvariableL <- length(between)
  withinvariableL <- length(within)
  varlabels <- toupper(colnames(data))
  if (!is.null(subjectid)) {
    data[,subjectid[1]] <- unlist(as.character(data[,subjectid[1]])) # force to characters
    indivparticipant <- sort(unique(as.character(data[,subjectid[1]])))
  } else {
    indivparticipant <- sort(unique(as.character(rownames(data))))
    if ((withinvariableL == 0) & (betweenvariableL > 0)) {
      # no subject ID parameter was provided and the analyis is between subjects
      # assume each is unique
      data$subjectid <- 1:nrow(data)
      subjectid <- "subjectid"
    } else {
      # a within subjects parameter has been specified
      stop("Error in RmimicAnova: A within subjects factor has been requested, but no subjectid parameter was provided to ensure pairwise comparisons. Please provide a column of subject ids.")
    }
  }
  
  completedata <- data[,c(subjectid[1], dependentvariable[1], between, within)]
  completedata <- completedata[complete.cases(completedata),] # remove missing datapoints
  
  # make sure data is collapsed across unnecessary observations
  tempcal <- sprintf("completedata <- doBy::summaryBy(%s ~", dependentvariable[1])
  faclist <- c(subjectid[1], between, within)
  for (cB in 1:length(faclist)) {
    if (cB == 1) {
      tempcal <- sprintf('%s %s', tempcal, faclist[cB])
    } else {
      tempcal <- sprintf('%s + %s', tempcal, faclist[cB])
    }
  }
  tempcal <- sprintf("%s, FUN=c(mean), data=completedata, keep.names=TRUE)", tempcal)
  suppressWarnings(eval(parse(text=tempcal)))
  rm(tempcal, cB)
  
  # Prepare output
  res <- list()
  
  # Output model
  demoout <- FALSE
  if ((verbose == TRUE) & (is.null(currentlevel))) {
    temptext <- "Univariate ANOVA Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    # Output model
    outstring <- sprintf("Analysis of %s were conducted using a", dependentvariable[1])
    # Number (Variable: Factors) 
    faclist <- c(between, within)
    for (cB in 1:length(faclist)) {
      factorname <- tolower(faclist[cB])
      factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
      factorsinvolved <- tolower(unique(unlist(as.character(completedata[,faclist[cB]]))))
      outstring <- sprintf("%s %d (%s: %s)", outstring, length(factorsinvolved), factorname, paste(factorsinvolved, collapse = ", "))
      if (cB < length(faclist)) {
        outstring <- sprintf("%s \u00D7", outstring) # seems to work on Mac as well
      }
    }
    outstring <- sprintf('%s univariate', outstring)
    if (!is.null(within)) {
      outstring <- sprintf('%s repeated measures', outstring)
    }
    
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s ANOVA using the ez (Lawrence, %s) and Rmimic (Pontifex, %s) packages in %s.', outstring, strsplit(as.character(packageDate("ez")),"-")[[1]][1],strsplit(as.character(packageDate("Rmimic")),"-")[[1]][1], rvers)
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    if (verbosedescriptives != FALSE) {
      demoout <- TRUE
    }
  } # end verbose
  
  # Compute demographics
  res$demographics <- Rmimic::descriptives(completedata[,c(dependentvariable[1], between, within)], groupvariable = c(between, within), verbose=demoout)
  
  # ezANOVA is not presently able to take inputs using dynamic variable names
  funcal <- sprintf('result <- ez::ezANOVA(data=completedata,dv=%s,wid=%s,', dependentvariable[1], subjectid[1])
  if (!is.null(between)) {
    tfuncal <- '.('
    for (cB in 1:betweenvariableL) {
      tfuncal <- sprintf('%s%s', tfuncal,between[cB])
      if (cB < betweenvariableL) {
        tfuncal <- sprintf('%s,', tfuncal)
      }
    }
    tfuncal <- sprintf('%s)', tfuncal)
    funcal <- sprintf('%sbetween=%s,', funcal, tfuncal)
    rm(tfuncal)
  } else {
    funcal <- sprintf('%sbetween=NULL,', funcal)
  }
  if (!is.null(within)) {
    tfuncal <- '.('
    for (cB in 1:withinvariableL) {
      tfuncal <- sprintf('%s%s', tfuncal,within[cB])
      if (cB < withinvariableL) {
        tfuncal <- sprintf('%s,', tfuncal)
      }
    }
    rm(cB)
    tfuncal <- sprintf('%s)', tfuncal)
    funcal <- sprintf('%swithin=%s,', funcal, tfuncal)
    rm(tfuncal)
  } else {
    funcal <- sprintf('%s, within=NULL,', funcal)
  }
  funcal <- sprintf('%stype=3,detailed=TRUE,return_aov=TRUE)', funcal)
  
  # Evaluate the text string and tell ezANOVA to shut up, outputs to results
  suppressWarnings(eval(parse(text=funcal)))
  #res$call <- funcal
  rm(funcal)
  
  # Trim off intercept from ANOVA and significance indicator
  result$ANOVA <- result$ANOVA[2:nrow(result$ANOVA),]
  #result$ANOVA <- result$ANOVA[,which(colnames(result$ANOVA) != "p<.05")]
  
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
      temp <- psychometric::CI.Rsq(tempR2, length(indivparticipant), length(unlist(strsplit(as.character(result$ANOVA$Effect[cT]),"[:]"))), level = confidenceinterval)
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
    temptext <- sprintf('%s %s) = ', temptext, pullvalue)
    # F stat
    temptext <- sprintf('%s%.1f, p', temptext, round(as.numeric(result$ANOVA$F[cR]), digits = 1))
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
  
  res$aov <- result$aov
  res$stats <- result$ANOVA
  
  if ((verbose == TRUE) & (is.null(currentlevel))) {
    
    # Show ANOVA Table
    if (!is.null(sphericity)) {
      summarylabel <- sprintf("Summary %s corrected ANOVA table", sphericity)
    } else {
      summarylabel <- "Summary ANOVA table"
    }
    cat(sprintf("\n"))
    Rmimic::typewriter(summarylabel, tabs=0, spaces=0, characters=floor(spansize*.9))
    if (toupper(feffect) == toupper("Partial Eta Squared")) {
      summarylabel <- sprintf("With partial eta squared based cohens f effect statistics.")
    } else {
      summarylabel <- sprintf("With generalized eta squared based cohens f effect statistics.")
    }
    Rmimic::typewriter(summarylabel, tabs=0, spaces=0, characters=floor(spansize*.9))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # create header labels
    vectnames <- c("Effect", "df", "F", "p")
    if (operatingsystem == "Windows") {
      temptext <- sprintf("f\u00b2 [%2.0f%% CI]", floor(confidenceinterval*100))
    } else {
      temptext <- sprintf("f^2 [%2.0f%% CI]", floor(confidenceinterval*100))
    }
    vectnames <- c(vectnames, temptext)
    sepgap <- matrix(floor(spansize/length(vectnames)), nrow=1, ncol=length(vectnames))
    sepgap[1,length(vectnames)-1] <- sepgap[1,length(vectnames)-1]-2
    sepgap[1,length(vectnames)] <- sepgap[1,length(vectnames)]-2
    colnames(sepgap) <- vectnames
    for (cC in 1:length(vectnames)) {
      if (vectnames[cC] == "Effect") {
        cat(sprintf("%-*s",sepgap[1,cC],vectnames[cC]))
      } else {
        cat(sprintf("%-*s",sepgap[1,cC],vectnames[cC]))
      }
    }
    cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    for (cR in 1:nrow(res$stats)) {
      for (cC in 1:length(vectnames)) {
        if (vectnames[cC] == "Effect") {
          cat(sprintf("%-*s\n",sepgap[1,cC],res$stats$Effect[cR]))
          cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:1]))), " "), collapse = "")))
        } else if (vectnames[cC] == "df") {
          pullval <- as.character(c(res$stats$DFn[cR], res$stats$DFd[cR]))
          pullval <- sprintf("%s, %s",pullval[1], pullval[2])
          cat(sprintf("%-*s",sepgap[1,cC],pullval))
        } else if (vectnames[cC] == "F") {
          cat(sprintf("%-*.1f",sepgap[1,cC],round(res$stats$F[cR],digits=1)))
        } else if (vectnames[cC] == "p") {
          # report P value
          outPvalue <- Rmimic::fuzzyP(res$stats$p[cR])
          if (outPvalue$modifier == "=") {
            pullvalue <- outPvalue$report
          } else {
            pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
          }
          if (outPvalue$interpret <= studywiseAlpha) {
            pullvalue <- sprintf('%s%s', pullvalue, "**")
          }
          cat(sprintf("%-*s",sepgap[1,cC],pullvalue))
          rm(outPvalue)
        } else {
          # effect size
          tempfsq <- sprintf("%.2f", round(as.numeric(res$stats$fsquared[cR]), digits = 2))
          if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
            tempfsq <- "< 0.01"
          }
          tempfsqlw <- sprintf("%.2f", round(as.numeric(res$stats$fsquared.ci.lower[cR]), digits = 2))
          if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
            tempfsqlw <- "0.0"
          }
          tempfsqup <- sprintf("%.2f", round(as.numeric(res$stats$fsquared.ci.upper[cR]), digits = 2))
          if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
            tempfsqup <- "0.0"
          }
          pullvalue <- sprintf('%s [%s, %s]', tempfsq, tempfsqlw, tempfsqup)
          cat(sprintf("%-*s",sepgap[1,cC], pullvalue))
        }
      } # end cC
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    } # end cR
  
    # output text writeups
    cat(sprintf("\n\n%s\n", "Summary ANOVA with posthoc decompositions, means, SD, and test statistics"))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  } # end verbose
  
  if (is.null(currentlevel)) {
    currentlevelout <- 0
  } else {
    currentlevelout <- currentlevel
  }
  
  # populate posthoc decompositions
  posthocANOVA <- NULL
  posthocttest <- NULL
  for (cR in 1:nrow(result$ANOVA)) {
    
    effectname <- result$ANOVA$Effect[cR]
    if ((verbose == TRUE) & (is.null(currentlevel))) {
      Rmimic::typewriter(effectname, tabs=currentlevelout, spaces=0, characters=floor(spansize*.9))
    }
    
    # check to see if planned contrast
    forcetrig <- 0
    if (!is.null(planned)) {
      if (planned == TRUE) {
        forcetrig <- 2
      } else {
        for (cXR in 1:length(planned)) {
          if (effectname == planned[cXR]) {
            forcetrig <- 1
          }
        }
        rm(cXR)
      }
    }
    
    # snag p value 
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$ANOVA$p[cR]))
    if (verbose == TRUE) {
      temp <- unlist(strsplit(as.character(effectname),"[:]"))
      outtext <- sprintf('There was')
      if (outPvalue$interpret <= studywiseAlpha) {
        if (length(temp) > 1) {
          outtext <- sprintf("%s an interaction of %s", outtext, paste(temp, collapse=sprintf(" \u00D7 ")))
        } else {
          outtext <- sprintf("%s a main effect of %s", outtext, temp[1])
        }
      } else {
        if (length(temp) > 1) {
          outtext <- sprintf("%s no interaction of %s", outtext, paste(temp, collapse=sprintf(" \u00D7 ")))
        } else {
          outtext <- sprintf("%s no main effect of %s", outtext, temp[1])
        }
      }
      outtext <- sprintf("%s, %s", outtext, result$ANOVA$textoutput[cR])
      Rmimic::typewriter(outtext, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
      rm(outtext)
      cat(sprintf("\n"))
    }
    
    if ((outPvalue$interpret <= studywiseAlpha) | (forcetrig > 0)) {
      # effect was significant or planned
      factorsinvolved <- unlist(strsplit(as.character(result$ANOVA$Effect[cR]),"[:]"))
      factorsinvolvedL <- length(factorsinvolved)
      if (factorsinvolvedL < 4) {
        
        # subset database for only those factors
        workingdatabase <- completedata[,c(subjectid[1], dependentvariable[1], factorsinvolved)]
        
        # make sure data is collapsed across unnecessary observations
        tempcal <- sprintf("workingdatabase <- doBy::summaryBy(%s ~", dependentvariable[1])
        faclist <- c(subjectid[1], factorsinvolved)
        for (cB in 1:length(faclist)) {
          if (cB == 1) {
            tempcal <- sprintf('%s %s', tempcal, faclist[cB])
          } else {
            tempcal <- sprintf('%s + %s', tempcal, faclist[cB])
          }
        }
        tempcal <- sprintf("%s, FUN=c(mean), data=workingdatabase, keep.names=TRUE)", tempcal)
        suppressWarnings(eval(parse(text=tempcal)))
        rm(tempcal, cB)
        
        subbetween <- NULL
        if (factorsinvolved %in% between) {
          subbetween <- factorsinvolved
        }
        subwithin <- NULL
        if (factorsinvolved %in% within) {
          subwithin <- factorsinvolved
        }
        
        if (factorsinvolvedL == 1) {
          # main effect
          
          # run t-test
          ttestresult <- Rmimic::RmimicTtest(workingdatabase, 
                dependentvariable=dependentvariable[1], 
                subjectid=subjectid[1], 
                between=subbetween, 
                within=subwithin,
                collapse=TRUE, nonparametric=nonparametric, posthoc=posthoc, 
                confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE)
          
          # modify output to indicate subtest
          for (cE in 1:nrow(ttestresult$stats)) {
            ttestresult$stats$Comparison[cE] <- sprintf('Main effect of %s; %s', as.character(result$ANOVA$Effect[cR]), ttestresult$stats$Comparison[cE])
            ttestresult$stats$EffectLabel <- as.character(result$ANOVA$Effect[cR])
          }
          
          # merge output 
          if (is.null(posthocttest)) {
            posthocttest <- ttestresult$stats
          } else {
            posthocttest <- rbind(posthocttest,ttestresult$stats)
          }
          
          # output findings
          if (verbose == TRUE) {
            # main effect decomposed with t-test
            Rmimic::typewriter(sprintf("%s Post hoc Comparisons %s", paste(replicate(5, spancharacter), collapse = ""), paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
            
            for (cE in 1:nrow(ttestresult$stats)) {
              Rmimic::typewriter(ttestresult$stats$interpretation[cE], tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              cat(sprintf("\n"))
            }
          }
          
          rm(subbetween,subwithin, ttestresult)
        } else {
          # interaction
          
          # see what we are working with
          factorlengthmatrix <- data.frame(matrix(NA,nrow=1, ncol=factorsinvolvedL))
          colnames(factorlengthmatrix) <- factorsinvolved
          for (cB in 1:factorsinvolvedL) {
            factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(workingdatabase[,factorsinvolved[cB]]))))
          }
          # deal with the most factors first
          factorlengthmatrix <- factorlengthmatrix[,order(-factorlengthmatrix[1,])]
          
          # loop through each factor
          for (cB in 1:factorsinvolvedL) {
            currentfactorinvolved <- colnames(factorlengthmatrix)[cB]
            currentfactorlevelsinvolved <- unique(unlist(as.character(workingdatabase[,currentfactorinvolved])))
            otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactorinvolved)]
            
            if (verbose == TRUE) {
              Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cB, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
            
              if (length(otherfactorsinvolved) == 1) {
                outval <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactorinvolved[1])
                Rmimic::typewriter(outval, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
                cat(sprintf('\n'))
              } else {
                outval <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactorinvolved[1])
                Rmimic::typewriter(outval, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
                cat(sprintf('\n'))
              }
            }
            
            # hold current factor level constant and subset
            for (cD in 1:length(currentfactorlevelsinvolved)) {
              subworkingdatabase <- workingdatabase[which(workingdatabase[,currentfactorinvolved] == currentfactorlevelsinvolved[cD]),]
              
              if (verbose == TRUE) {
                outval <- sprintf('For %s: %s', currentfactorinvolved, currentfactorlevelsinvolved[cD])
                Rmimic::typewriter(outval, tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              }
              
              # make sure data is collapsed across unnecessary observations
              tempcal <- sprintf("subworkingdatabase <- doBy::summaryBy(%s ~", dependentvariable[1])
              faclist <- c(subjectid[1], otherfactorsinvolved)
              for (cE in 1:length(faclist)) {
                if (cE == 1) {
                  tempcal <- sprintf('%s %s', tempcal, faclist[cE])
                } else {
                  tempcal <- sprintf('%s + %s', tempcal, faclist[cE])
                }
              }
              tempcal <- sprintf("%s, FUN=c(mean), data=subworkingdatabase, keep.names=TRUE)", tempcal)
              suppressWarnings(eval(parse(text=tempcal)))
              rm(tempcal, cE)
              
              subbetween <- NULL
              subwithin <- NULL
              for (cE in 1:length(otherfactorsinvolved)) {
                if (otherfactorsinvolved[cE] %in% between) {
                  subbetween <- otherfactorsinvolved[cE]
                }
                if (otherfactorsinvolved[cE] %in% within) {
                  subwithin <- otherfactorsinvolved[cE]
                }
              }
              rm(cE)
              
              booltrig <- 0
              # more than one additional variable
              if (length(otherfactorsinvolved) > 1) { 
                booltrig <- 1
              } else {
                # one factor but with more than 2 levels
                otherfactorlevelsinvolved <- unique(unlist(as.character(subworkingdatabase[,otherfactorsinvolved[1]])))
                if (length(otherfactorlevelsinvolved) > 2) { 
                  booltrig <- 1
                }
              }
              
              if (booltrig == 0) {
                # only a single factor with less than 3 levels
                # only a t-test is required
                
                # run t-test
                ttestresult <- Rmimic::RmimicTtest(subworkingdatabase, 
                       dependentvariable=dependentvariable[1], 
                       subjectid=subjectid[1], 
                       between=subbetween, 
                       within=subwithin,
                       collapse=TRUE, nonparametric=nonparametric, posthoc=posthoc, 
                       confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE)

                # modify output to indicate subtest
                for (cE in 1:nrow(ttestresult$stats)) {
                  ttestresult$stats$Comparison[cE] <- sprintf('Interaction of %s; %s: %s, %s', as.character(result$ANOVA$Effect[cR]), currentfactorinvolved, currentfactorlevelsinvolved[cD], ttestresult$stats$Comparison[cE])
                  ttestresult$stats$EffectLabel <- as.character(result$ANOVA$Effect[cR])
                }
                
                # merge output 
                if (is.null(posthocttest)) {
                  posthocttest <- ttestresult$stats
                } else {
                  posthocttest <- rbind(posthocttest,ttestresult$stats)
                }
                
                # output findings
                if (verbose == TRUE) {
                  for (cE in 1:nrow(ttestresult$stats)) {
                    Rmimic::typewriter(ttestresult$stats$interpretation[cE], tabs=currentlevelout+3, spaces=0, characters=floor(spansize*.9))
                    cat(sprintf("\n"))
                  }
                }
                
                rm(subbetween,subwithin, ttestresult, subworkingdatabase)
                
              } else {
                # need to run full anova again
                
                posthocplanned <- NULL
                if (forcetrig > 0) {
                  posthocplanned <- TRUE
                }
                subresult <- RmimicAnova(subworkingdatabase,
                       dependentvariable=dependentvariable[1], 
                       subjectid=subjectid[1], 
                       between=subbetween, 
                       within=subwithin,
                       currentlevel=(currentlevelout+2),
                       nonparametric=nonparametric, posthoc=posthoc,
                       sphericity=sphericity, planned=posthocplanned,
                       verbose=TRUE)
                
                # modify ANOVA output to indicate subtest
                for (cE in 1:nrow(subresult$stats)) {
                  subresult$stats$Effect[cE] <- sprintf('Interaction of %s; %s: %s, %s', as.character(result$ANOVA$Effect[cR]), currentfactorinvolved, currentfactorlevelsinvolved[cD], subresult$stats$Effect[cE])
                  subresult$stats$EffectLabel <- as.character(result$ANOVA$Effect[cR])
                }
                if (is.null(posthocANOVA)) {
                  posthocANOVA <- subresult$stats
                } else {
                  posthocANOVA <- rbind(posthocANOVA,subresult$stats)
                }
                
                # modify Posthoc ttest output to indicate subtest
                for (cE in 1:nrow(subresult$posthocttest)) {
                  subresult$posthocttest$Comparison[cE] <- sprintf('Interaction of %s; %s: %s, %s', as.character(result$ANOVA$Effect[cR]), currentfactorinvolved, currentfactorlevelsinvolved[cD], subresult$posthocttest$Comparison[cE])
                  subresult$posthocttest$EffectLabel <- as.character(result$ANOVA$Effect[cR])
                }
                if (is.null(posthocttest)) {
                  posthocttest <- subresult$posthocttest
                } else {
                  posthocttest <- rbind(posthocttest,subresult$posthocttest)
                }
                
                # output findings
                if (verbose == TRUE) {
                  
                }
                cat(sprintf("\n"))
              }
              rm (subworkingdatabase, booltrig)
            }
            rm(cD)
          }
        }
      } else {
        # more than 3 factors involved
        
        outtext <- sprintf("An interaction exceeding 3 variables was detected. To conserve computational resouces, interactions exceeding 3 variables should be decomposed in a stepwise fashion manually.")
        Rmimic::typewriter(outtext, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
        rm(outtext)
        
        # loop through each factor
        for (cB in 1:factorsinvolvedL) {
          currentfactorinvolved <- factorsinvolved[cB]
          otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactorinvolved)]
          
          Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cB, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
          
          outval <- sprintf("Post-hoc decomposition of the %s interaction should be conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactorinvolved[1])
          Rmimic::typewriter(outval, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
          cat(sprintf("\n"))
        } # end cB
        
      } # less than 4 factors involved
    } # effect is significant
  } # for each effect
  res$posthocANOVA <- posthocANOVA
  res$posthocttest <- posthocttest
  
  if ((verbose == TRUE) & (is.null(currentlevel))) {
    cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  options(warn = oldw) # turn warnings back to original settings
  return(res)
}
  