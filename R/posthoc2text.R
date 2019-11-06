#' posthoc2text
#'
#' @description Output ANOVA posthoc results in APA style format with effect sizes and confidence intervals.
#'
#' @param result Anova and posthoc test output
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param spansize Number of characters to include on a line.
#' @param currentlevelout Integer for tab indenting used in recursion.
#' @param posthoclimit Parameter to specify the limit for breaking down interaction terms. Default is 6 indicating a 6 way interaction would not be automatically broken down.
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 30, 2019
#' 
#' @export

posthoc2text <- function(result, studywiseAlpha=0.05, spansize=95, currentlevelout=0, posthoclimit=6) {
  
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
    
    # if posthoc tests were run
    if (outPvalue$interpret <= studywiseAlpha) {
      if (result$stats$EffectPostHoc[cR] == 1) {
        if (length(temp) == 1) {
          # must be a main effect
          # look for associated t tests
          if ("posthocttest" %in% names(result)) {
            # subset posthoc tests for this effect
            subposthocttest <- result$posthocttest[which(result$posthocttest$EffectNumber == cR),]
            if (nrow(subposthocttest) > 0) {
              numswitchtext <- "Comparison"
              if (nrow(subposthocttest) > 1) {
                numswitchtext <- sprintf('%ss', numswitchtext)
              }
              Rmimic::typewriter(sprintf("%s Post hoc %s %s", paste(replicate(5, spancharacter), collapse = ""), numswitchtext, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              for (cE in 1:nrow(subposthocttest)) {
                Rmimic::typewriter(subposthocttest$interpretation[cE], tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
                cat(sprintf("\n"))
              }
            }
            rm(subposthocttest, cE)
          }
        } else {
          # must be an interaction
          
          if (length(temp) > posthoclimit) {
            
            outtext <- sprintf("For clarity purposes, interactions exceeding %d variables are not automatically printed. The Rmimic::posthoc2text() function can be used to show these results.\n", (posthoclimit-1))
            Rmimic::typewriter(outtext, tabs=currentlevelout+1, spaces=0, characters=floor(spansize*.9))
            rm(outtext)
            
          } else {
            
            # reconstruct ANOVA levels
            factorsinvolved <- unlist(strsplit(as.character(result$stats$Effect[cR]),"[:]"))
            factorsinvolvedL <- length(factorsinvolved)
            factorlengthmatrix <- data.frame(matrix(NA,nrow=2, ncol=factorsinvolvedL))
            colnames(factorlengthmatrix) <- factorsinvolved
            effectlevels <- trimws(unlist(strsplit(substr(result$stats$EffectLevels[cR],2,nchar(result$stats$EffectLevels[cR])-1),"[:]")))
            for (currentLevels in 1:length(effectlevels)) {
              tempvect <- trimws(unlist(strsplit(effectlevels[currentLevels], "[(]")))
              factorlengthmatrix[1,trimws(unlist(strsplit(tempvect[1], " ")))[1]] <- trimws(unlist(strsplit(tempvect[1], " ")))[2]
              factorlengthmatrix[2,trimws(unlist(strsplit(tempvect[1], " ")))[1]] <- substr(tempvect[2],1,nchar(tempvect[2])-1)
            }
            
            # loop through each factor involved
            for (cB in 1:factorsinvolvedL) {
              currentfactorinvolved <- factorsinvolved[cB]
              currentfactorlevelsinvolved <- trimws(unlist(strsplit(factorlengthmatrix[2, currentfactorinvolved], ",")))
              otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactorinvolved)]
              
              if (length(otherfactorsinvolved) == 1) {
                decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactorinvolved[1])
              } else {
                decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactorinvolved[1])
              }
              Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cB, paste(replicate(5, spancharacter), collapse = "")), tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              Rmimic::typewriter(decomptext, tabs=currentlevelout+2, spaces=0, characters=floor(spansize*.9))
              cat(sprintf("\n"))
              decompdir <- paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 "))
              
              # hold current factor level constant and subset
              for (cD in 1:length(currentfactorlevelsinvolved)) {
                decompconst <- sprintf('For %s: %s', currentfactorinvolved, currentfactorlevelsinvolved[cD])
                Rmimic::typewriter(decompconst, tabs=currentlevelout+3, spaces=0, characters=floor(spansize*.9))
                
                # determine what needed to be done for this breakdown
                booltrig <- 0
                # more than one additional variable
                if (length(otherfactorsinvolved) > 1) { 
                  booltrig <- 2
                } else {
                  # one factor but with more than 2 levels
                  if (as.numeric(factorlengthmatrix[1,otherfactorsinvolved]) > 2) { 
                    booltrig <- 1
                  }
                }
                
                if (booltrig == 0) {
                  # only a single factor with less than 3 levels
                  # only a t-test is required
                
                  # subset posthoc tests for this effect
                  if ("posthocttest" %in% names(result)) {
                    subposthocttest <- result$posthocttest[which(result$posthocttest$EffectDirection == sprintf('[%s] %s', decompconst, decompdir)),]
                    if (nrow(subposthocttest) > 0) {
                      for (cE in 1:nrow(subposthocttest)) {
                        Rmimic::typewriter(subposthocttest$interpretation[cE], tabs=currentlevelout+4, spaces=0, characters=floor(spansize*.9))
                        cat(sprintf("\n"))
                      }
                    }
                    rm(subposthocttest, cE)
                  }
                } else {
                  # a full ANVOA should have been run
                  # these now exist as seperate entries
                  
                  outlabel <- sprintf('PosthocANOVA_%s_%s_%s', paste(factorsinvolved, collapse="By"), currentfactorlevelsinvolved[cD], paste(otherfactorsinvolved, collapse="By"))
                  if (outlabel %in% names(result)) {
                  
                    funcal <- sprintf('backflipanovaout <- result$%s', outlabel)
                    suppressWarnings(eval(parse(text=funcal)))
                    
                    # perform backflip
                    posthoc2text(backflipanovaout, studywiseAlpha=studywiseAlpha, spansize=spansize, currentlevelout=currentlevelout+3)
                    
                  }
                } # end t test or anova
              } # end cD
            } # end cB
          }
        } # end interaction check
      } # end posthoc check
    } # was significant
    cat(sprintf("\n"))
  } # end each row of stats
}
  
  
  