#' posthoc2text
#'
#' @description Output ANOVA posthoc results in APA style format with effect sizes and confidence intervals.
#'
#' @param result Anova and posthoc test output
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param spansize Number of characters to include on a line.
#' @param currentlevelout Integer for tab indenting used in recursion.
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 30, 2019
#' 
#' @export

posthoc2text <- function(result, studywiseAlpha=0.05, spansize=95, currentlevelout=0) {
  
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

  for (cR in 1:nrow(result$stats)) {
    effectname <- result$stats$Effect[cR]
    if (currentlevelout == 0) {
      cat(sprintf('%s\n', effectname))
    }
    
    temp <- unlist(strsplit(as.character(effectname),"[:]"))
    outtext <- sprintf('There was')
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$stats$p[cR]))
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
    outtext <- sprintf("%s, %s", outtext, result$stats$textoutput[cR])
    Rmimic::typewriter(outtext, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
    rm(outtext)
    cat(sprintf("\n"))
    
    if (length(temp) == 1) {
      # must be a main effect
      # look for associated t tests
      if ("posthocttest" %in% names(result)) {
        # subset posthoc tests for this effect
        subposthocttest <- result$posthocttest[which(result$posthocttest$EffectLabel == effectname),]
        if (nrow(subposthocttest) > 0) {
          Rmimic::typewriter(sprintf("%s Post hoc Comparisons %s", paste(replicate(5, spancharacter), collapse = ""), paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
          for (cE in 1:nrow(subposthocttest)) {
            Rmimic::typewriter(subposthocttest$interpretation[cE], tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
            cat(sprintf("\n"))
          }
        }
        rm(subposthocttest, cE)
      }
    } else {
      # must be an interaction
      
      # loop through each factor involved
      for (cB in 1:length(temp)) {
        currentfactorinvolved <- temp[cB]
        otherfactorsinvolved <- temp[which(temp != currentfactorinvolved)]
        decompdir <- c(currentfactorinvolved, otherfactorsinvolved)
        decompdir <- paste(decompdir, collapse=sprintf(" \u00D7 "))
        
        # determine what was done first
        boolANOVA <- FALSE
        hitlistANOVA <- c()
        boolTTEST <- FALSE
        hitlistTTEST <- c()
        if ("posthocANOVA" %in% names(result)) {
          for (cA in 1:nrow(result$posthocANOVA)) {
            tempcounteffectlabel <- unlist(strsplit(as.character(result$posthocANOVA$EffectLabel[cA]),"[;]"))
            tempcounteffectdirection <- unlist(strsplit(as.character(result$posthocANOVA$EffectDirection[cA]),"[;]"))
            if ((tempcounteffectlabel[1] == effectname) & (tempcounteffectdirection[1] == decompdir)) {
              # matches first pull
              boolANOVA <- TRUE
              hitlistANOVA <- c(hitlistANOVA, cA)
            }
          }
        }
        if ("posthocttest" %in% names(result)) {
          for (cA in 1:nrow(result$posthocttest)) {
            tempcounteffectlabel <- unlist(strsplit(as.character(result$posthocttest$EffectLabel[cA]),"[;]"))
            tempcounteffectdirection <- unlist(strsplit(as.character(result$posthocttest$EffectDirection[cA]),"[;]"))
            if ((tempcounteffectlabel[1] == effectname) & (tempcounteffectdirection[1] == decompdir)) {
              # matches first pull
              boolTTEST <- TRUE
              hitlistTTEST <- c(hitlistTTEST, cA)
            }
          }
        }
        
        if ((boolTTEST == TRUE) & (boolANOVA == FALSE)) {
          # Interaction must have been able to be broken down with just t-tests
          Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cB, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
          subposthocttest <- result$posthocttest[hitlistTTEST,]
          
          # check for multiple approaches
          approachlist <- unique(unlist(as.character(subposthocttest$Decomptext)))
          for (cA in 1:length(approachlist)) {
            Rmimic::typewriter(approachlist[cA], tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
            cat(sprintf('\n'))
            workingsubposthocttest <- subposthocttest[which(subposthocttest$Decomptext == approachlist[cA]),]
            
            for (cAS in 1:nrow(workingsubposthocttest)) {
              outval <- sprintf('For %s: %s', currentfactorinvolved, workingsubposthocttest$DecompFor[cAS])
              Rmimic::typewriter(outval, tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              Rmimic::typewriter(workingsubposthocttest$interpretation[cAS], tabs=currentlevelout+3, spaces=0, characters=floor(spansize*.9))
              cat(sprintf('\n'))
            }
          }
        } else if (boolANOVA == TRUE) {
          # Interaction must have required another anova
          Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cB, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
          subposthocANOVA <- result$posthocANOVA[hitlistANOVA,]
          if (boolTTEST == TRUE) {
            # one of the ANOVA effects must have been broken down
            subposthocttest <- result$posthocttest[hitlistTTEST,]
          }
          # check for multiple approaches
          approachlist <- c()
          for (cA in 1:nrow(subposthocANOVA)) {
            tempval <- unlist(strsplit(as.character(subposthocANOVA$Decomptext[cA]),"[;]"))
            approachlist <- c(approachlist, tempval[1])
          }
          for (cA in 1:length(unique(approachlist))) {
            Rmimic::typewriter(unique(approachlist)[cA], tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
            cat(sprintf('\n'))
            workingsubposthocANOVA <- subposthocANOVA[which(approachlist == unique(approachlist)[cA]),]
            forlist <- c()
            for (cF in 1:nrow(workingsubposthocANOVA)) {
              tempval <- unlist(strsplit(as.character(workingsubposthocANOVA$DecompFor[cF]),"[;]"))
              forlist <- c(forlist, tempval[1])
            }
            for (cF in 1:length(unique(forlist))) {
              outval <- sprintf('For %s: %s', currentfactorinvolved, unique(forlist)[cF])
              Rmimic::typewriter(outval, tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              forworkingsubposthocANOVA <- workingsubposthocANOVA[which(forlist == unique(forlist)[cF]),]
              
              # strip off preceeding breakdowns
              for (cFA in 1:nrow(forworkingsubposthocANOVA)) {
                workinglabels <- c('EffectLabel', 'EffectDirection', 'DecompFor', 'Decomptext')
                for (cFAL in 1:length(workinglabels)) {
                  tempval <- unlist(strsplit(as.character(forworkingsubposthocANOVA[cFA, workinglabels[cFAL]]),"[;]"))
                  if (length(tempval) == 1) {
                    tempvaltext <- forworkingsubposthocANOVA$Effect[cFA]
                  } else if (length(tempval) == 2) {
                    tempvaltext <- trimws(tempval[2])
                  } else if (length(tempval) > 2) {
                    tempvaltext <- trimws(paste(tempval[2:length(tempval)], collapse=";"))
                  }
                  if (length(tempvaltext) > 0) {
                    forworkingsubposthocANOVA[cFA, workinglabels[cFAL]] <- tempvaltext
                  } else {
                    forworkingsubposthocANOVA[cFA, workinglabels[cFAL]] <- NA
                  }
                  
                }
              }
              backflipanovaout <- list()
              backflipanovaout$stats <- forworkingsubposthocANOVA
              rm(forworkingsubposthocANOVA)
              
              if (boolTTEST == TRUE) {
                hitlistTTEST <- c()
                for (cFA in 1:nrow(subposthocttest)) {
                  tempcounteffectlabel <- unlist(strsplit(as.character(subposthocttest$Decomptext[cFA]),"[;]"))
                  tempcounteffectdirection <- unlist(strsplit(as.character(subposthocttest$DecompFor[cFA]),"[;]"))
                  if ((tempcounteffectlabel[1] == unique(approachlist)[cA]) & (tempcounteffectdirection[1] == unique(forlist)[cF])) {
                    hitlistTTEST <- c(hitlistTTEST, cFA)
                  }
                }
                forsubposthocttest <- subposthocttest[hitlistTTEST,]
                
                # strip off preceeding breakdowns
                for (cFA in 1:nrow(forsubposthocttest)) {
                  workinglabels <- c('EffectLabel', 'EffectDirection', 'DecompFor', 'Decomptext')
                  for (cFAL in 1:length(workinglabels)) {
                    tempval <- unlist(strsplit(as.character(forsubposthocttest[cFA, workinglabels[cFAL]]),"[;]"))
                    if (length(tempval) == 1) {
                      tempvaltext <- forsubposthocttest$Effect[cFA]
                    } else if (length(tempval) == 2) {
                      tempvaltext <- trimws(tempval[2])
                    } else if (length(tempval) > 2) {
                      tempvaltext <- trimws(paste(tempval[2:length(tempval)], collapse=";"))
                    }
                    if (length(tempvaltext) > 0) {
                      forsubposthocttest[cFA, workinglabels[cFAL]] <- tempvaltext
                    } else {
                      forsubposthocttest[cFA, workinglabels[cFAL]] <- NA
                    }
                  }
                }
                backflipanovaout$posthocttest <- forsubposthocttest
                rm(forsubposthocttest)
              }
              
              # this is where you could do a backflip into this same function!
              posthoc2text(backflipanovaout, studywiseAlpha=studywiseAlpha, spansize=spansize, currentlevelout=currentlevelout+2)
              
            }
          }
        }
      }
    } # end main effect vs interaction
    cat(sprintf("\n"))
  } # end each row of stats
}
  
  
  