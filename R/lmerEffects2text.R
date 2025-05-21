#' lmerEffects2text
#'
#' @description Sub function to clean up and round values from Rmimic::lmerEffects.
#'
#' @param res list containing data frames for stats, randomstats, and rsquared values
#' @param tag string containing ID tag for the text output
#'
#' @return res list with an added column for the text output
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, March 28, 2025
#' 
#' @importFrom stringr str_replace_all str_remove_all str_split
#'
#' @export

lmerEffects2text <- function(res, tag='', subtag='', testconfidence=FALSE, significanceconfidence=FALSE) {
  
  studywiseAlpha <- tryCatch({
    studywiseAlpha <- res$studywiseAlpha
  }, error = function(e) {
    studywiseAlpha <- NULL
  })
  if (is.null(studywiseAlpha)) {
    studywiseAlpha <- 0.05
  }
  
  confidenceinterval <- tryCatch({
    confidenceinterval <- res$confidenceinterval
  }, error = function(e) {
    confidenceinterval <- 0.95
  })
  
  # clean up values for reporting
  tempdbs <- tryCatch({
    tempdbs <- res$stats
  }, error = function(e) {
    tempdbs <- NULL
  })
  if (!is.null(tempdbs)) {
    dataframeout <- data.frame(matrix(NA,nrow=nrow(tempdbs),ncol=11))
    colnames(dataframeout) <- c('Effect', 'df', 'F', 'p', 'pinterp', 'fsquared', 'fsquaredCI', 'textoutput', 'idtag', 'significance', 'factorsinvolved')
    dataframeout$Effect <- tempdbs$Effect
    
    for (cR in 1:nrow(tempdbs)) {
      # DF
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$DFn[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      temptext <- sprintf('%s', pullvalue)
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$DFd[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      temptext <- sprintf('%s, %s', temptext, pullvalue)
      dataframeout$df[cR] <- temptext
      
      #f
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$F.value[cR]), digits = 1))
      if (pullvalue == "0.0") {
        pullvalue = "&lt; 0.1"
      }
      dataframeout$F[cR] <- pullvalue
      
      # P val
      outPvalue <- Rmimic::fuzzyP(as.numeric(tempdbs$p.value[cR]), studywiseAlpha=studywiseAlpha, html=TRUE)
      temptext <- sprintf('%s', outPvalue$report)
      if (outPvalue$modifier != '&#61;') {
        temptext <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      dataframeout$p[cR] <- tempdbs$p.value[cR]
      dataframeout$pinterp[cR] <- outPvalue$interpret
      dataframeout$significance[cR] <- outPvalue$significance
      
      # textoutput
      temptext <- sprintf('<span style="font-style: italic;">F</span>(%s)', dataframeout$df[cR])
      pullvalue <- dataframeout$F[cR]
      if (substr(pullvalue, 1, 1) == "&lt;") {
        temptext <- sprintf('%s %s', temptext, pullvalue)
      } else {
        temptext <- sprintf('%s &#61; %s', temptext, pullvalue)
      }
      
      if (testconfidence) {
        if (!is.na(tempdbs$F.value.ci.lower[cR])) {
          tempfsqlw <- sprintf("%.1f", round(as.numeric(tempdbs$F.value.ci.lower[cR]), digits = 1))
          if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
            tempfsqlw <- "0.0"
          }
          tempfsqup <- sprintf("%.1f", round(as.numeric(tempdbs$F.value.ci.upper[cR]), digits = 1))
          if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
            tempfsqup <- "0.0"
          }
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
      
      temptext <- sprintf('%s, <span style="font-style: italic;">p</span>', temptext)
      pullvalue <- Rmimic::fuzzyP(dataframeout$p[cR], studywiseAlpha=studywiseAlpha, html=TRUE)
      if ((pullvalue$report == "0.000") | (pullvalue$report == "0.00") | (pullvalue$report == "0.0") | (pullvalue$report == "0") | (pullvalue$report == "0.")) {
        pullvalue$report <- "0.001"
        pullvalue$modifier <- "&lt;"
      }
      temptext <- sprintf('%s %s %s', temptext, pullvalue$modifier, pullvalue$report)
      
      if (significanceconfidence) {
        if (!is.na(tempdbs$p.value.ci.lower[cR])) {
          tempfsqlw <- sprintf("%.3f", round(as.numeric(tempdbs$p.value.ci.lower[cR]), digits = 3))
          if ((tempfsqlw == "-0.000") | (tempfsqlw == "0.000") | (tempfsqlw == "-0.00") | (tempfsqlw == "0.00") | (tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
            #tempfsqlw <- "0.0"
            #if (pullvalue$modifier == '&lt;') {
              tempfsqlw <- "&lt; 0.001"
            #}
          }
          tempfsqup <- sprintf("%.3f", round(as.numeric(tempdbs$p.value.ci.upper[cR]), digits = 3))
          if ((tempfsqup == "-0.000") | (tempfsqup == "0.000") | (tempfsqup == "-0.00") | (tempfsqup == "0.00") | (tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
            #tempfsqup <- "0.0"
            #if (pullvalue$modifier == '&lt;') {
              tempfsqup <- "&lt; 0.001"
            #}
          }
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
      
      # effect size
      tempfsq <- sprintf("%.2f", round(as.numeric(tempdbs$fsquared[cR]), digits = 2))
      if ((tempfsq == "-0.00") | (tempfsq == "0.00")) {
        tempfsq <- "&lt; 0.01"
      }
      dataframeout$fsquared[cR] <- tempfsq
      
      # effect size ci
      tempfsqlw <- sprintf("%.2f", round(as.numeric(tempdbs$fsquared.ci.lower[cR]), digits = 2))
      if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
        tempfsqlw <- "0.0"
      }
      tempfsqup <- sprintf("%.2f", round(as.numeric(tempdbs$fsquared.ci.upper[cR]), digits = 2))
      if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
        tempfsqup <- "0.0"
      }
      if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
        # Negative
        if (outPvalue$significance) {
          if (tempfsqup == "0.0") {
            # cant really be significant and zero - likely rounding issue
            tempfsqup <- "-0.01"
          }
        }
      } else {
        # Positive
        if (outPvalue$significance) {
          if (tempfsqlw == "0.0") {
            # cant really be significant and zero - likely rounding issue
            tempfsqlw <- "0.01"
          }
        }
      }
      dataframeout$fsquaredCI[cR] <- sprintf('%s, %s', tempfsqlw, tempfsqup)
      
      temptext <- sprintf('%s, <span id="effectsizestatistic"><span style="font-style: italic;">f\u200a\u00b2</span> &#61; %s [%2.0f%% CI: %s to %s]</span>.', temptext,
                          dataframeout$fsquared[cR], floor(confidenceinterval*100),
                          tempfsqlw, tempfsqup)
      dataframeout$textoutput[cR] <- temptext
      
      
      factorsinvolved <- unlist(strsplit(as.character(dataframeout$Effect[cR]),"[:]"))
      factorsinvolvedL <- length(factorsinvolved)
      if (outPvalue$significance) {
        if (factorsinvolvedL == 1) {
          textout <- 'There was a main effect of'
        } else {
          textout <- 'There was an interaction of'
        }
      } else {
        if (factorsinvolvedL == 1) {
          textout <- 'There was no main effect of'
        } else {
          textout <- 'There was no interaction between'
        }
      }
      grouptext <- paste0(factorsinvolved, sep="", collapse = " &times; ")
      dataframeout$factorsinvolved[cR] <- factorsinvolvedL
      
      effecttext <- stringr::str_replace_all(as.character(dataframeout$Effect[cR]), '[:]', '_by_')
      if (tag == '') {
        idtextout <- sprintf('fixed-%s', effecttext)
      } else {
        idtextout <- sprintf('%s-fixed-%s', tag, effecttext)
      }
      
      dataframeout$textoutput[cR] <- sprintf('<span id="%s">%s %s, %s</span>', idtextout, textout, grouptext, temptext)
      dataframeout$idtag[cR] <- idtextout 
    }
    res$stats$textoutput <- dataframeout$textoutput
    res$stats$idtag <- dataframeout$idtag
    res$stats$significance <- dataframeout$significance
    res$stats$factorsinvolved <- dataframeout$factorsinvolved
  }
  
  
  # clean up values for reporting
  tempdbs <- tryCatch({
    tempdbs <- res$randomstats
  }, error = function(e) {
    tempdbs <- NULL
  })
  if (!is.null(tempdbs)) {
    dataframeout <- data.frame(matrix(NA,nrow=nrow(tempdbs),ncol=9))
    colnames(dataframeout) <- c('Effect', 'df', 'Likelihood Ratio', 'log-likelihood', 'p', 'pinterp', 'textoutput', 'idtag', 'significance')
    dataframeout$Effect <- tempdbs$Effect
    
    for (cR in 1:nrow(tempdbs)) {
      # DF
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$DF[cR]), digits = 1))
      if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
        if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
          pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
        }
      }
      dataframeout$df[cR] <- sprintf('%s', pullvalue)
      
      #Likelihood Ratio
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$LRT[cR]), digits = 1))
      if (pullvalue == "0.0") {
        pullvalue = "&lt; 0.1"
      }
      dataframeout$'Likelihood Ratio'[cR] <- pullvalue
      
      #log-likelihood
      pullvalue <- sprintf('%.1f', round(as.numeric(tempdbs$LogLikelihood[cR]), digits = 1))
      if (pullvalue == "0.0") {
        pullvalue = "&lt; 0.1"
      }
      dataframeout$'log-likelihood'[cR] <- pullvalue
      
      # P val
      outPvalue <- Rmimic::fuzzyP(as.numeric(tempdbs$p.value[cR]), studywiseAlpha=studywiseAlpha, html=TRUE)
      dataframeout$pinterp[cR] <- outPvalue$interpret
      dataframeout$significance[cR] <- outPvalue$significance
      
      # text out
      temptext <- sprintf('<span style="font-style: italic;">Likelihood Ratio</span> (%s)', dataframeout$df[cR])
      pullvalue <- dataframeout$`Likelihood Ratio`[cR]
      if (substr(pullvalue, 1, 1) == "&lt;") {
        temptext <- sprintf('%s %s, ', temptext, pullvalue)
      } else {
        temptext <- sprintf('%s &#61; %s, ', temptext, pullvalue)
      }
      temptext <- sprintf('%s<span style="font-style: italic;">log-likelihood</span>', temptext)
      pullvalue <- dataframeout$`log-likelihood`[cR]
      if (substr(pullvalue, 1, 1) == "&lt;") {
        temptext <- sprintf('%s %s, ', temptext, pullvalue)
      } else {
        temptext <- sprintf('%s &#61; %s, ', temptext, pullvalue)
      }
      temptext <- sprintf('%s<span style="font-style: italic;">p</span>', temptext)
      if ((outPvalue$report == "0.000") | (outPvalue$report == "0.00") | (outPvalue$report == "0.0") | (outPvalue$report == "0") | (outPvalue$report == "0.")) {
        outPvalue$report <- "0.001"
        outPvalue$modifier <- "&lt;"
      }
      temptext <- sprintf('%s %s %s.', temptext, outPvalue$modifier, outPvalue$report)
      
      
      effecttext <- stringr::str_remove_all(as.character(dataframeout$Effect[cR]), '[(]')
      effecttext <- stringr::str_remove_all(effecttext, '[)]')
      effecttext <- stringr::str_replace_all(effecttext, stringr::fixed(" "), "")
      effecttext <- stringr::str_replace_all(effecttext, '[|]', '_within_')
      
      if (tag == '') {
        idtextout <- sprintf('random-%s', effecttext)
      } else {
        idtextout <- sprintf('%s-random-%s', tag, effecttext)
      }
      
      dataframeout$textoutput[cR] <- sprintf('<span id="%s">%s</span>', idtextout, temptext)
      dataframeout$idtag[cR] <- idtextout 
      
    }
    res$randomstats$textoutput <- dataframeout$textoutput
    res$randomstats$idtag <- dataframeout$idtag
    res$randomstats$significance <- dataframeout$significance
  }
  
  
  # clean up values for reporting
  tempdbs <- tryCatch({
    tempdbs <- res$rsquared
  }, error = function(e) {
    tempdbs <- NULL
  })
  if (!is.null(tempdbs)) {
    dataframeout <- data.frame(matrix(NA,nrow=nrow(tempdbs),ncol=5))
    colnames(dataframeout) <- c('Portion', 'Effect', 'CI', 'textoutput', 'idtag')
    dataframeout$Portion <- tempdbs$portion
    
    for (cR in 1:nrow(tempdbs)) {
      #effect
      pullvalue <- sprintf('%.2f', round(as.numeric(tempdbs$effects[cR]), digits = 2))
      if ((pullvalue == "-0.00") | (pullvalue == "0.00")) {
        pullvalue = "&lt; 0.1"
      } else {
        pullvalue = sprintf("&#61; %s", pullvalue)
      }
      dataframeout$Effect[cR] <- pullvalue
      
      # effect size ci
      tempfsqlw <- sprintf("%.2f", round(as.numeric(tempdbs$ci.lower[cR]), digits = 2))
      if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
        tempfsqlw <- "0.0"
      }
      tempfsqup <- sprintf("%.2f", round(as.numeric(tempdbs$ci.upper[cR]), digits = 2))
      if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
        tempfsqup <- "0.0"
      }
      dataframeout$CI[cR] <- sprintf('%s to %s', tempfsqlw, tempfsqup)
      
      temptext <- sprintf('<span style="font-style: italic;">R\u200a\u00b2</span> %s [%2.0f%% CI: %s]</span>.',
                          dataframeout$Effect[cR],
                          floor(confidenceinterval*100),
                          dataframeout$CI[cR])
      
      if (tag == '') {
        idtextout <- sprintf('rsquared-%s', effecttext)
      } else {
        idtextout <- sprintf('%s-rsquared-%s', tag, effecttext)
      }
      
      dataframeout$textoutput[cR] <- sprintf('<span id="%s">%s</span>', idtextout, temptext)
      dataframeout$idtag[cR] <- idtextout 
      
    }
    
    res$rsquared$textoutput <- dataframeout$textoutput
    res$rsquared$idtag <- dataframeout$idtag
  }
  
  
  tempdbs <- tryCatch({
    tempdbs <- res$posthoc
  }, error = function(e) {
    tempdbs <- NULL
  })
  if (!is.null(tempdbs)) {
  
    outputnames <- names(tempdbs)
    for (cOutputNames in 1:length(outputnames)) {
      workingdataout <- NULL
      textcall <- sprintf("workingdataout <- tempdbs$%s", outputnames[cOutputNames])
      eval(parse(text=textcall))
      if (!is.null(workingdataout)) {
        if (is.data.frame(workingdataout)) {
          # posthoc tests
          
          workingdataout$textoutput <-NA
          for (cR in 1:nrow(workingdataout)) {
          
            temptext <- ''
            groupdescript <- sprintf('%s', workingdataout$C1name[cR])
            if (subtag != '') {
              groupdescript <- sprintf('%s<span class="explainmethodtext"><sub>%s</sub></span>', groupdescript, subtag)
            }
            groupdescript <- sprintf('%s (%.1f &plusmn; %.1f)', groupdescript, workingdataout$C1mean[cR], workingdataout$C1sd[cR])
            groupdescript <- sprintf('%s and %s', groupdescript, workingdataout$C2name[cR])
            if (subtag != '') {
              groupdescript <- sprintf('%s<span class="explainmethodtext"><sub>%s</sub></span>', groupdescript, subtag)
            }
            groupdescript <- sprintf('%s (%.1f &plusmn; %.1f)', groupdescript, workingdataout$C2mean[cR], workingdataout$C2sd[cR])
            
            if (workingdataout$significant[cR]) {
              temptext <- sprintf('The difference between %s was statistically significant;', groupdescript)
            } else {
              temptext <- sprintf('No significant differences were observed between %s;', groupdescript)
            }
          
            pullvalue <- sprintf('%.1f', workingdataout$t.ratio[cR])
            if ((pullvalue == "0.0") | (pullvalue == "-0.0")) {
              pullvalue = "&lt; 0.1"
            } else {
              pullvalue <- sprintf('&#61; %s', pullvalue)
            }
            if (is.infinite(workingdataout$df[cR])) {
              temptext <- sprintf('%s <span style="font-style: italic;">t</span>(&infin;) %s', temptext, pullvalue)
            } else {
              temptext <- sprintf('%s <span style="font-style: italic;">t</span>(%d) %s', temptext, floor(workingdataout$df[cR]), pullvalue)
            }
            
            if (testconfidence) {
              if (!is.na(workingdataout$t.conf.int.lower[cR])) {
                tempfsqlw <- sprintf("%.1f", round(as.numeric(workingdataout$t.conf.int.lower[cR]), digits = 1))
                if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
                  tempfsqlw <- "0.0"
                }
                tempfsqup <- sprintf("%.1f", round(as.numeric(workingdataout$t.conf.int.upper[cR]), digits = 1))
                if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
                  tempfsqup <- "0.0"
                }
                temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                    floor(confidenceinterval*100), tempfsqlw, tempfsqup)
              }
            }
            
            outPvalue <- Rmimic::fuzzyP(as.numeric(workingdataout$p.value[cR]), studywiseAlpha=studywiseAlpha, html=TRUE)
            if ((outPvalue$report == "0.000") | (outPvalue$report == "0.00") | (outPvalue$report == "0.0") | (outPvalue$report == "0") | (outPvalue$report == "0.")) {
              outPvalue$report <- "0.001"
              outPvalue$modifier <- "&lt;"
            }
            pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
            temptext <- sprintf('%s, <span style="font-style: italic;">p</span> %s', temptext, pullvalue)
            
            if (significanceconfidence) {
              if (!is.na(workingdataout$p.conf.int.lower[cR])) {
                tempfsqlw <- sprintf("%.3f", round(as.numeric(workingdataout$p.conf.int.lower[cR]), digits = 3))
                if ((tempfsqlw == "-0.000") | (tempfsqlw == "0.000")) {
                  tempfsqlw <- "0.0"
                  if (outPvalue$modifier == '&lt;') {
                    tempfsqlw <- "&lt; 0.001"
                  }
                }
                tempfsqup <- sprintf("%.3f", round(as.numeric(workingdataout$p.conf.int.upper[cR]), digits = 3))
                if ((tempfsqup == "-0.000") | (tempfsqup == "0.000")) {
                  tempfsqup <- "0.0"
                  if (outPvalue$modifier == '&lt;') {
                    tempfsqup <- "&lt; 0.001"
                  }
                }
                temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                    floor(confidenceinterval*100), tempfsqlw, tempfsqup)
              }
            }
            
            # effect size
            if (!is.na(workingdataout$correlation[cR])) {
              # within subjects
              effectsizetext <- sprintf('<span style="font-style: italic;">d<sub class="sub">rm</sub></span>')
            } else{
              # between subjects
              effectsizetext <- sprintf('<span style="font-style: italic;">d<sub class="sub">s</sub></span>')
            }
            
            tempfsq <- sprintf("%.2f", workingdataout$effectsize[cR])
            if ((tempfsq == "-0.00") | (tempfsq == "0.00")) {
              tempfsq <- "&lt; 0.01"
            } else {
              tempfsq <- sprintf("&#61; %s", tempfsq)
            }
            
            # effect size ci
            tempfsqlw <- sprintf("%.2f", workingdataout$effectsize.conf.int.lower[cR])
            if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
              tempfsqlw <- "0.0"
            }
            tempfsqup <- sprintf("%.2f", workingdataout$effectsize.conf.int.upper[cR])
            if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
              tempfsqup <- "0.0"
            }
            if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
              # Negative
              if (workingdataout$significant[cR]) {
                if (tempfsqup == "0.0") {
                  # cant really be significant and zero - likely rounding issue
                  tempfsqup <- "-0.01"
                }
              }
            } else {
              # Positive
              if (workingdataout$significant[cR]) {
                if (tempfsqlw == "0.0") {
                  # cant really be significant and zero - likely rounding issue
                  tempfsqlw <- "0.01"
                }
              }
            }
            
            effectsizetext <- sprintf('%s %s [%2.0f%% CI: %s to %s]', effectsizetext, tempfsq, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
            temptext <- sprintf('%s, <span id=\"effectsizestatistic\">%s</span>.', temptext, effectsizetext)
            
            workingdataout$textoutput[cR] <- sprintf('<span>%s</span>', temptext)
          } # for each row of workingdataout
          
          # put data back in
          textcall <- sprintf("res$posthoc$%s$textoutput <- workingdataout$textoutput", outputnames[cOutputNames])
          eval(parse(text=textcall))
          
          # see if we need to get fancy
          if (nrow(workingdataout) == 1) {
            # only a single posthoc test was required
            
            # get effect localizer
            tempname <- stringr::str_split(outputnames[cOutputNames], 'Posthoc_')[[1]][2]
            tempname <- stringr::str_split(tempname, '_for_')[[1]][1]
            ftesteffecttext <- res$stats$textoutput[which(res$stats$Effect == tempname)]
            
            # swap cohens f for cohens d given only two levels
            ftesteffecttext <- stringr::str_split(ftesteffecttext, '<span id=\"effectsizestatistic\">')[[1]]
            ttesteffecttext <- stringr::str_split(workingdataout$textoutput[1], '<span id=\"effectsizestatistic\">')[[1]]
            ftesteffecttext <- sprintf('%s<span id=\"effectsizestatistic\">%s', ftesteffecttext[1], ttesteffecttext[2])
            
            # put effect back into the table
            res$stats$textoutput[which(res$stats$Effect == tempname)] <- ftesteffecttext
          }
          
        } else {
          # another ANOVA
          
          # reprocess recursively
          workingdataout <- lmerEffects2text(workingdataout, tag=tag, subtag=subtag, testconfidence=testconfidence, significanceconfidence=significanceconfidence)
          
          # put it back
          textcall <- sprintf("res$posthoc$%s <- workingdataout", outputnames[cOutputNames])
          eval(parse(text=textcall))
        }
      }
    }
  }
  
  return(res)
}

