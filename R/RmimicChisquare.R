#' RmimicChisquare
#'
#' @description Compute SPSS style results for Chi-square analysis with effect size and confidence intervals. For samples less than 1000, Fishers exact test statistic is used.
#'
#' @param data Data frame containing the variables of interest.
#' @param variables Variable name or list of variables to compute descriptives for.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param planned Boolean operator for if the posthoc tests should be outputted regardless of the test statistic. Default is FALSE.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#'
#' @return
#' \item{stats}{Summary of analysis.}
#' \item{tabularoutput}{Chi square table.}
#' \item{table}{Simple frequency table.}
#' \item{posthocttest}{A summary table for posthoc tests.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 6, 2020
#'
#' @importFrom epitools oddsratio.fisher oddsratio.wald
#'
#' @export

RmimicChisquare <- function(variables=FALSE, data=FALSE, confidenceinterval=0.95, studywiseAlpha=0.05, planned=FALSE, verbose=TRUE) {  
  
  # debug
  #variables <- c('MotherBMI', 'ChildBMI')
  #data <- workingdata
  #variables=c('carb','cyl')
  #data=mtcars
  #variables=c('treatment','improvement')
  #data=workingdata
  #confidenceinterval=0.95
  #studywiseAlpha=0.05
  #verbose=TRUE
  #planned=TRUE

  
  # Prepare output
  res <- list()
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

  if (variables[1] == FALSE) {
    variables <- names(data)
    cDF <- data
  } else {
    cDF <- data.frame(data[,variables])
    names(cDF) <- variables
  }
  cDF <- cDF[which(complete.cases(cDF)),]
  
  # make sure things are factors
  if (!is.factor(cDF[,1])) {
    cDF[,1] <- factor(cDF[,1])
  }
  if (!is.factor(cDF[,1])) {
    cDF[,2] <- factor(cDF[,2])
  }
  
  # kick to subprocess
  res <- Rmimic::subprocessRmimicChisquare(variables=variables, data=cDF)
  
  boolposthoc <- TRUE
  results <- list()
  results$data <- FALSE
  if (nrow(cDF) < 1000) {
    tryCatch(results <- epitools::oddsratio.fisher(x=cDF[,1], y=cDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval), error=function(e){boolposthoc <- FALSE})
  } else {
    tryCatch(results <- epitools::oddsratio.wald(x=cDF[,1], y=cDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval), error=function(e){boolposthoc <- FALSE})
  }
  if (results$data[1] == FALSE) {
    boolposthoc <- FALSE
  }
  
  if (boolposthoc == TRUE) {
    res$table <- results$data # place data table in output
  
    # populate posthoc decompositions
    predictorlevels <- levels(cDF[,1])
    predictorlevelsL <- length(predictorlevels)
    outcomelevels <- levels(cDF[,2])
    outcomelevelsL <- length(outcomelevels)
    
    masterdataframeout <- data.frame(matrix(NA,nrow=0,ncol=10))
    colnames(masterdataframeout) <- c('Reference','Predictor','Outcome','OddsRatio', 'OR.lower.conf.int','OR.upper.conf.int', 'p.value', 'textoutput','overallmodel.report', 'overallmodel.p.value')
    
    # outcomes should binary to make the most sense
    for (cOutcomeLevelsB in 1:(outcomelevelsL)) {
      for (cOutcomeLevelsC in 1:(outcomelevelsL)) {
        if ((cOutcomeLevelsB != cOutcomeLevelsC) & (cOutcomeLevelsB/cOutcomeLevelsC <= 1)) {
      
          # create temporary dataset containing only the binary outcome contrast
          posthoccDF1 <- cDF[which(cDF[,2] == outcomelevels[cOutcomeLevelsB]),]
          posthoccDF2 <- cDF[which(cDF[,2] == outcomelevels[cOutcomeLevelsC]),]
          posthoccDF <- rbind(posthoccDF1,posthoccDF2)
          rm(posthoccDF1,posthoccDF2)
          # rerun on subset
          tempout <- subprocessRmimicChisquare(variables=names(posthoccDF), data=posthoccDF)
          
          # each predictor gets a turn being the reference
          for (cPredictorLevelsB in 1:(predictorlevelsL)) {
            
            dataframeout <- data.frame(matrix(NA,nrow=(predictorlevelsL-1),ncol=10))
            colnames(dataframeout) <- c('Outcome','Reference','Predictor','OddsRatio', 'OR.lower.conf.int','OR.upper.conf.int', 'p.value', 'textoutput', 'overallmodel.report', 'overallmodel.p.value')
            
            subposthoccDF <- posthoccDF
            baselevel <- predictorlevels[cPredictorLevelsB]
            subposthoccDF[,1] <- factor(subposthoccDF[,1], levels = c(baselevel, predictorlevels[which(predictorlevels != baselevel)])) 
            
            if (nrow(subposthoccDF) < 1000) {
              results <- epitools::oddsratio.fisher(x=subposthoccDF[,1], y=subposthoccDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval)
            } else {
              results <- epitools::oddsratio.wald(x=subposthoccDF[,1], y=subposthoccDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval)
            }
            
            # cycle through output
            for (cR in 2:(nrow(results$measure))) {
              dataframeout$Outcome[cR-1] <- sprintf('%s vs %s', outcomelevels[cOutcomeLevelsB], outcomelevels[cOutcomeLevelsC])
              dataframeout$Reference[cR-1] <- baselevel
              dataframeout$Predictor[cR-1] <- rownames(results$measure)[cR]
              dataframeout$OddsRatio[cR-1] <- results$measure[cR,1]
              dataframeout$OR.lower.conf.int[cR-1] <- results$measure[cR,2]
              dataframeout$OR.upper.conf.int[cR-1] <- results$measure[cR,3]
              if (nrow(subposthoccDF) < 1000) {
                dataframeout$p.value[cR-1] <- results$p.value[cR,2]
              } else {
                dataframeout$p.value[cR-1] <- results$p.value[cR,3]
              }
              outtextstring <- sprintf("Odds Ratio = %.2f", round(dataframeout$OddsRatio[cR-1], digits=2))
              outtextstring <- sprintf("%s [%2.0f%% CI: %.2f to %.2f]", outtextstring, confidenceinterval*100, dataframeout$OR.lower.conf.int[cR-1], dataframeout$OR.upper.conf.int[cR-1])
              outPvalue <- Rmimic::fuzzyP(as.numeric(dataframeout$p.value[cR-1]))
              outtextstring <- sprintf('%s, p %s %s', outtextstring, outPvalue$modifier, outPvalue$report)
              dataframeout$textoutput[cR-1] <- outtextstring
              dataframeout$overallmodel.report[cR-1] <- tempout$stats$textoutput[1]
              dataframeout$overallmodel.p.value[cR-1] <- tempout$stats$p.value[1]
            }
            
            masterdataframeout <- rbind(masterdataframeout, dataframeout)
          
          } # end loop cPredictorLevelsB
        } # prevent duplicates
      } # end loop cOutcomeLevelsC
    } # end loop cOutcomeLevelsB
    res$posthocttest <- masterdataframeout
    
  } # posthoc tests got borked
    
  if (verbose == TRUE) {
    
    temptext <- "Chi-square Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- "Analysis were conducted using the"
    outstring <- sprintf('%s gmodels (Warnes, Bolker, Lumley, & Johnson, %s)', outstring, strsplit(as.character(utils::packageDate("gmodels")),"-")[[1]][1])
    outstring <- sprintf('%s, epitools (Aragon, %s)', outstring, strsplit(as.character(utils::packageDate("epitools")),"-")[[1]][1])
    outstring <- sprintf('%s, and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    if (nrow(cDF) < 1000) {
      outstring <- sprintf('%s Test statistics and odds ratios were computed using %s Exact test and conditional maximum likelihood estimation.', outstring, sprintf("%s's", 'Fisher'))
    }
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # show model
    for (cR in 1:length(res$tabularoutput)) {
      cat(sprintf("%s\n",res$tabularoutput[cR]))
    }
    
    outstring <- "Interpretation with Test Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # report main outcome
    outtext <- ""
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$stats$p.value[1]))
    if (outPvalue$interpret <= studywiseAlpha) {
      outtext <- sprintf("The likelihood of %s was observed to relate to %s", variables[2], variables[1])
    } else {
      outtext <- sprintf("There was no difference in the expected frequency of %s as a function of %s", variables[2], variables[1])
    }
    outtext <- sprintf("%s, %s.", outtext, res$stats$textoutput[1])
    Rmimic::typewriter(outtext, tabs=0, spaces=2, characters=floor(spansize*.9))
    rm(outtext)
    cat(sprintf("\n"))
    
    if (boolposthoc == TRUE) {
    
      if ((outPvalue$interpret <= studywiseAlpha) | (planned == TRUE)) {
        
        tablevel <- 1
        posthocoutcomelevels <- unique(res$posthocttest$Outcome)
        if (length(posthocoutcomelevels) > 1) {
          tablevel <- 2
        }
        for (cPosthocoutcomelevels in 1:length(posthocoutcomelevels)) {
          subposthocDF <- res$posthocttest
          subposthocDF <- subposthocDF[which(subposthocDF$Outcome == posthocoutcomelevels[cPosthocoutcomelevels]),]
          tempvect <- strsplit(posthocoutcomelevels[cPosthocoutcomelevels], " ")[[1]]
          booltrig <- 1
          if (length(posthocoutcomelevels) > 1) {
            booltrig <- 0
            
            Rmimic::typewriter(sprintf('Outcome %s: %s', variables[2], posthocoutcomelevels[cPosthocoutcomelevels]), tabs=1, spaces=0, characters=floor(spansize*.9))
            Rmimic::typewriter(paste(replicate(35, spancharacter), collapse = ""), tabs=1, spaces=0, characters=floor(spansize*.9))
            outtext <- "For clarity, analysis were repeated within each combination of paired outcome levels.\n"
            Rmimic::typewriter(outtext, tabs=1, spaces=0, characters=floor(spansize*.9))
            
            outPvalue2 <- Rmimic::fuzzyP(as.numeric(subposthocDF$overallmodel.p.value[1]))
            if (outPvalue2$interpret <= studywiseAlpha) {
              booltrig <- 1
              outtext <- sprintf("The likelihood of %s (%s, %s) was observed to relate to %s", variables[2], tempvect[1], tempvect[3], variables[1])
            } else {
              outtext <- sprintf("There was no difference in the expected frequency of %s (%s, %s) as function of %s", variables[2], tempvect[1], tempvect[3], variables[1])
            }
            outtext <- sprintf("%s, %s.", outtext, subposthocDF$overallmodel.report[1])
            Rmimic::typewriter(outtext, tabs=1, spaces=2, characters=floor(spansize*.9))
            rm(outtext)
            cat(sprintf("\n"))
          }
          
          if ((booltrig == 1) | (planned == TRUE)) {
              
            posthocpredictorlevels <- unique(subposthocDF$Reference)
            for (cPosthocpredictorlevels in 1:length(posthocpredictorlevels)) {
            
              subsubposthocDF <- subposthocDF[which(subposthocDF$Reference == posthocpredictorlevels[cPosthocpredictorlevels]),]
              
              Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cPosthocpredictorlevels, paste(replicate(5, spancharacter), collapse = "")), tabs=tablevel, spaces=2, characters=floor(spansize*.9))
              decomptext <- sprintf("Post-hoc comparisons were conducted by examining the expected frequency of %s (%s, %s) relative to %s: %s.\n", variables[2], tempvect[1], tempvect[3], variables[1], posthocpredictorlevels[cPosthocpredictorlevels])
              Rmimic::typewriter(decomptext, tabs=tablevel, spaces=2, characters=floor(spansize*.9))
              
              for (cR in 1:nrow(subsubposthocDF)) {
                
                directionalstate <- "an increased"
                if (subsubposthocDF$OddsRatio[cR] <= 1.0) {
                  directionalstate <- "a decreased"
                }
                outPvalue3 <- Rmimic::fuzzyP(as.numeric(subsubposthocDF$p.value[cR]))
                
                tempinsttext <- sprintf('%s %s',posthocpredictorlevels[cPosthocpredictorlevels], variables[1])
                outtext <- sprintf("%s %s", subsubposthocDF$Predictor[cR], variables[1]) 
                if (outPvalue3$interpret <= studywiseAlpha) {
                  outtext <- sprintf("%s was associated with %s likelihood, relative to %s, of %s %s",outtext, directionalstate,tempinsttext,tempvect[3], variables[2])
                } else {
                  outtext <- sprintf("%s made no difference, relative to %s, in the expected frequency of %s %s",outtext,tempinsttext, tempvect[3], variables[2])
                }
                outtext <- sprintf("%s, %s.", outtext, subsubposthocDF$textoutput[cR])
                Rmimic::typewriter(outtext, tabs=tablevel+1, spaces=2, characters=floor(spansize*.9))
                rm(outtext)
                cat(sprintf("\n"))
                
              }
              
              cat(sprintf("\n"))
            } # end loop cPosthocpredictorlevels
          } # booltrig or planned
        } # end loop cPosthocoutcomelevels
      } # if significant or planned
    } # if posthoc test were run
    
    # end output
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  } # end verbose
  
  #return(res)
}
  
  