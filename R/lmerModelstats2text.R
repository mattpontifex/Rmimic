#' lmerModelstats2text
#'
#' @description Sub function to create text interepretations from Rmimic::lmerModelstats
#'
#' @param res list containing output from Rmimic::lmerModelstats
#'
#' @return res list with an added elements for the text output
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, June 3, 2025
#'
#' @export

lmerModelstats2text <- function(res, confidenceinterval=0.95, studywiseAlpha=0.05, testconfidence=FALSE, significanceconfidence=FALSE) {
  
  temptext <- ""
  
  if ('Chisq' %in% names(res)) {
    if ('df' %in% names(res)) {
      temptext <- sprintf('%s<span style="font-style: italic;">&chi;&hairsp;\u00b2</span>', temptext)
      temptext <- sprintf('%s (%d)', temptext, ceiling(res$df))
      temptext <- sprintf('%s = %.1f', temptext, recursiveround(res$Chisq, digits=1))
      
      if (testconfidence) {
        if ('Chisq.ci.lower' %in% names(res)) {
          tempfsqlw <- sprintf("%.1f", recursiveround(res$Chisq.ci.lower, digits=1))
          if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
            tempfsqlw <- "0.0"
          }
          tempfsqup <- sprintf("%.1f", recursiveround(res$Chisq.ci.upper, digits=1))
          if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
            tempfsqup <- "0.0"
          }
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
    }
  }
  
  if ('p.value' %in% names(res)) {
    pullvalue <- Rmimic::fuzzyP(res$p.value, studywiseAlpha=studywiseAlpha, html=TRUE)
    if (nchar(temptext) > 0) {
      temptext <- sprintf('%s, ', temptext)
    }
    temptext <- sprintf('%s<span style="font-style: italic;">p</span>', temptext)
    temptext <- sprintf('%s %s %s', temptext, pullvalue$modifier, pullvalue$report)
    
    if (significanceconfidence) {
      if ('p.value.ci.lower' %in% names(res)) {
        tempfsqlw <- sprintf("%.3f", recursiveround(res$p.value.ci.lower, digits=3))
        tempfsqup <- sprintf("%.3f", recursiveround(res$p.value.ci.upper, digits=3))
        temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                            floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
    }
  }
  
  if ('fsquared' %in% names(res)) {
    if (nchar(temptext) > 0) {
      temptext <- sprintf('%s, ', temptext)
    }
    temptext <- sprintf('%s<span style="font-style: italic;">f&hairsp;&hairsp;\u00b2</span>', temptext)
    temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$fsquared, digits=2))
    tempfsqlw <- sprintf("%.2f", recursiveround(res$fsquared.ci.lower, digits=2))
    if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
      tempfsqlw <- "0.0"
    }
    tempfsqup <- sprintf("%.2f", recursiveround(res$fsquared.ci.upper, digits=2))
    if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
      tempfsqup <- "0.0"
    }
    if ('p.value' %in% names(res)) {
      if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
        # Negative
        if (pullvalue$significance) {
          if (tempfsqup == "0.0") {
            # cant really be significant and zero - likely rounding issue
            tempfsqup <- "-0.01"
          }
        }
      } else {
        # Positive
        if (pullvalue$significance) {
          if (tempfsqlw == "0.0") {
            # cant really be significant and zero - likely rounding issue
            tempfsqlw <- "0.01"
          }
        }
      }
    }
    temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
  }
  
  if ('fixed.r.squared' %in% names(res)) {
    if (nchar(temptext) > 0) {
      temptext <- sprintf('%s, ', temptext)
    }
    temptext <- sprintf('%s<span style="font-style: italic;">R&hairsp;\u00b2<sub>fixed</sub></span>', temptext)
    temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$fixed.r.squared, digits=2))
    
    if (testconfidence) {
      if ('fixed.r.squared.ci.lower' %in% names(res)) {
        tempfsqlw <- sprintf("%.2f", recursiveround(res$fixed.r.squared.ci.lower, digits=2))
        if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
          tempfsqlw <- "0.0"
        }
        tempfsqup <- sprintf("%.2f", recursiveround(res$fixed.r.squared.ci.upper, digits=2))
        if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
          tempfsqup <- "0.0"
        }
        temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                            floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
    }
    
    temptext <- sprintf('%s, <span style="font-style: italic;">R&hairsp;\u00b2<sub>random</sub></span>', temptext)
    temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$random.r.squared, digits=2))
    
    if (testconfidence) {
      if ('random.r.squared.ci.lower' %in% names(res)) {
        tempfsqlw <- sprintf("%.2f", recursiveround(res$random.r.squared.ci.lower, digits=2))
        if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
          tempfsqlw <- "0.0"
        }
        tempfsqup <- sprintf("%.2f", recursiveround(res$random.r.squared.ci.upper, digits=2))
        if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
          tempfsqup <- "0.0"
        }
        temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                            floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
    }
    
  }
  
  if ('VIF' %in% names(res)) {
    if (nchar(temptext) > 0) {
      temptext <- sprintf('%s, ', temptext)
    }
    temptext <- sprintf('%sVIF = %.1f', temptext, recursiveround(res$VIF, digits=1))
    
    if (res$VIF > 2.5) {
      temptext <- sprintf("%s\n<br><sub>[Warning: VIF suggests a", temptext)
      if (res$VIF > 2.5) {
        degree <- "slight"
      }
      if (res$VIF > 5) {
        degree <- "moderate"
      }
      if (res$VIF > 10) {
        degree <- "severe"
      }
      temptext <- sprintf("%s %s", temptext, degree)
      temptext <- sprintf("%s multicollinearity issue]</sub>", temptext)
    }
  }
  
  res$text <- temptext
  #htmlTable::htmlTable(temptext)
  
  # manage change
  if ('change' %in% names(res)) {
    #res <- lmerModelstats(fit, altfit, show=FALSE)
    
    temptext <- ""
    if ('Chisq' %in% names(res$change)) {
      if ('df' %in% names(res$change)) {
        temptext <- sprintf('%s<span style="font-style: italic;">&chi;&hairsp;\u00b2<sub>change</sub></span>', temptext)
        temptext <- sprintf('%s (%d)', temptext, ceiling(res$change$df))
        temptext <- sprintf('%s = %.1f', temptext, recursiveround(res$change$Chisq, digits=1))
        
        if (testconfidence) {
          if ('Chisq.ci.lower' %in% names(res$change)) {
            tempfsqlw <- sprintf("%.1f", recursiveround(res$change$Chisq.ci.lower, digits=1))
            if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
              tempfsqlw <- "0.0"
            }
            tempfsqup <- sprintf("%.1f", recursiveround(res$change$Chisq.ci.upper, digits=1))
            if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
              tempfsqup <- "0.0"
            }
            temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                floor(confidenceinterval*100), tempfsqlw, tempfsqup)
          }
        }
        
      }
    }
    
    if ('p.value' %in% names(res$change)) {
      pullvalue <- Rmimic::fuzzyP(res$change$p.value, studywiseAlpha=studywiseAlpha, html=TRUE)
      if (nchar(temptext) > 0) {
        temptext <- sprintf('%s, ', temptext)
      }
      temptext <- sprintf('%s<span style="font-style: italic;">p</span>', temptext)
      temptext <- sprintf('%s %s %s', temptext, pullvalue$modifier, pullvalue$report)
      
      if (significanceconfidence) {
        if ('p.value.ci.lower' %in% names(res$change)) {
          tempfsqlw <- sprintf("%.3f", recursiveround(res$change$p.value.ci.lower, digits=3))
          tempfsqup <- sprintf("%.3f", recursiveround(res$change$p.value.ci.upper, digits=3))
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
    }
    
    if ('fsquared' %in% names(res$change)) {
      if (nchar(temptext) > 0) {
        temptext <- sprintf('%s, ', temptext)
      }
      temptext <- sprintf('%s<span style="font-style: italic;">f&hairsp;&hairsp;\u00b2<sub>change</sub></span>', temptext)
      temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$change$fsquared, digits=2))
      tempfsqlw <- sprintf("%.2f", recursiveround(res$change$fsquared.ci.lower, digits=2))
      if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
        tempfsqlw <- "0.0"
      }
      tempfsqup <- sprintf("%.2f", recursiveround(res$change$fsquared.ci.upper, digits=2))
      if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
        tempfsqup <- "0.0"
      }
      if ('p.value' %in% names(res$change)) {
        if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
          # Negative
          if (pullvalue$significance) {
            if (tempfsqup == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqup <- "-0.01"
            }
          }
        } else {
          # Positive
          if (pullvalue$significance) {
            if (tempfsqlw == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqlw <- "0.01"
            }
          }
        }
      }
      temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
    }
    
    if ('fixed.r.squared' %in% names(res$change)) {
      if (nchar(temptext) > 0) {
        temptext <- sprintf('%s, ', temptext)
      }
      temptext <- sprintf('%s<span style="font-style: italic;">&Delta;R&hairsp;\u00b2<sub>fixed</sub></span>', temptext)
      temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$change$fixed.r.squared, digits=2))
      
      if (testconfidence) {
        if ('fixed.r.squared.ci.lower' %in% names(res$change)) {
          tempfsqlw <- sprintf("%.2f", recursiveround(res$change$fixed.r.squared.ci.lower, digits=2))
          if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
            tempfsqlw <- "0.0"
          }
          tempfsqup <- sprintf("%.2f", recursiveround(res$change$fixed.r.squared.ci.upper, digits=2))
          if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
            tempfsqup <- "0.0"
          }
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
      
      temptext <- sprintf('%s, <span style="font-style: italic;">&Delta;R&hairsp;\u00b2<sub>random</sub></span>', temptext)
      temptext <- sprintf('%s = %.2f', temptext, recursiveround(res$change$random.r.squared, digits=2))
      
      if (testconfidence) {
        if ('random.r.squared.ci.lower' %in% names(res$change)) {
          tempfsqlw <- sprintf("%.2f", recursiveround(res$change$random.r.squared.ci.lower, digits=2))
          if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
            tempfsqlw <- "0.0"
          }
          tempfsqup <- sprintf("%.2f", recursiveround(res$change$random.r.squared.ci.upper, digits=2))
          if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
            tempfsqup <- "0.0"
          }
          temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                              floor(confidenceinterval*100), tempfsqlw, tempfsqup)
        }
      }
    }
    res$changetext <- temptext
  }
  
  # manage coefficients
  if ('coefficients' %in% names(res)) {
    res$coefficients$text <- NA
    for (cR in 1:nrow(res$coefficients)) {
      temptext <- ''
      #temptext <- sprintf("%s%s", temptext, res$coefficients$Variable[cR])
      temptext <- sprintf('%s<span style="font-style: italic;">B</span> = %.2f', temptext, recursiveround(res$coefficients$B[cR], digits=2))
      
      if ((!is.na(res$coefficients$B.lower.conf.int[cR])) & (!is.na(res$coefficients$B.upper.conf.int[cR]))) {
        tempfsqlw <- sprintf("%.2f", recursiveround(res$coefficients$B.lower.conf.int[cR], digits=2))
        if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
          tempfsqlw <- "0.0"
        }
        tempfsqup <- sprintf("%.2f", recursiveround(res$coefficients$B.upper.conf.int[cR], digits=2))
        if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
          tempfsqup <- "0.0"
        }
        pullvalue <- Rmimic::fuzzyP(res$coefficients$p.value[cR], studywiseAlpha=studywiseAlpha, html=TRUE)
        if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
          # Negative
          if (pullvalue$significance) {
            if (tempfsqup == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqup <- "-0.01"
            }
          }
        } else {
          # Positive
          if (pullvalue$significance) {
            if (tempfsqlw == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqlw <- "0.01"
            }
          }
        }
        temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
      temptext <- sprintf('%s, <span style="font-style: italic;">SE B</span> = %.2f', temptext, recursiveround(res$coefficients$SE[cR], digits=2))
      
      
      if (!is.na(res$coefficients$OddsRatio[cR])) {
        temptext <- sprintf('%s, <span style="font-style: italic;">Odds Ratio</span> = %.2f', temptext, recursiveround(res$coefficients$OddsRatio[cR], digits=2))
      } else {
        temptext <- sprintf('%s, <span style="font-style: italic;">&beta;</span> = %.2f', temptext, recursiveround(res$coefficients$Beta[cR], digits=2))
      }
      tempfsqlw <- NA
      tempfsqup <- NA
      if (!is.na(res$coefficients$OddsRatio[cR])) {
        if (!is.na(res$coefficients$OddsRatio.lower.conf.int[cR])) {
          tempfsqlw <- sprintf("%.2f", recursiveround(res$coefficients$OddsRatio.lower.conf.int[cR], digits=2))
        }
        if (!is.na(res$coefficients$OddsRatio.upper.conf.int[cR])) {
          tempfsqup <- sprintf("%.2f", recursiveround(res$coefficients$OddsRatio.upper.conf.int[cR], digits=2))
        }
      } else {
        if (!is.na(res$coefficients$Beta.lower.conf.int[cR])) {
          tempfsqlw <- sprintf("%.2f", recursiveround(res$coefficients$Beta.lower.conf.int[cR], digits=2))
        }
        if (!is.na(res$coefficients$Beta.upper.conf.int[cR])) {
          tempfsqup <- sprintf("%.2f", recursiveround(res$coefficients$Beta.upper.conf.int[cR], digits=2))
        }
      }
      if ((!is.na(tempfsqlw)) & (!is.na(tempfsqup))) {
        if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
          tempfsqlw <- "0.0"
        }
        if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
          tempfsqup <- "0.0"
        }
        pullvalue <- Rmimic::fuzzyP(res$coefficients$p.value[cR], studywiseAlpha=studywiseAlpha, html=TRUE)
        if (abs(as.numeric(tempfsqlw)) > abs(as.numeric(tempfsqup))) {
          # Negative
          if (pullvalue$significance) {
            if (tempfsqup == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqup <- "-0.01"
            }
          }
        } else {
          # Positive
          if (pullvalue$significance) {
            if (tempfsqlw == "0.0") {
              # cant really be significant and zero - likely rounding issue
              tempfsqlw <- "0.01"
            }
          }
        }
        temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext, floor(confidenceinterval*100), tempfsqlw, tempfsqup)
      }
      
      if (!is.na(res$coefficients$t[cR])) {
        temptext <- sprintf('%s, <span style="font-style: italic;">t</span>&hairsp;(%d) = %.1f', temptext, ceiling(res$coefficients$df[cR]),
                          recursiveround(res$coefficients$t[cR], digits=1))
      
        if (testconfidence) {
          if ('t.ci.lower' %in% names(res$coefficients)) {
            if ((!is.na(res$coefficients$t.ci.lower[cR])) & (!is.na(res$coefficients$t.ci.upper[cR]))) {
              tempfsqlw <- sprintf("%.1f", recursiveround(res$coefficients$t.ci.lower[cR], digits=1))
              if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
                tempfsqlw <- "0.0"
              }
              tempfsqup <- sprintf("%.1f", recursiveround(res$coefficients$t.ci.upper[cR], digits=1))
              if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
                tempfsqup <- "0.0"
              }
              temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                  floor(confidenceinterval*100), tempfsqlw, tempfsqup)
            }
          }
        }
      }
      if (!is.na(res$coefficients$z[cR])) {
        temptext <- sprintf('%s, <span style="font-style: italic;">z</span> = %.1f', temptext, 
                            recursiveround(res$coefficients$z[cR], digits=1))
        
        if (testconfidence) {
          if ('z.ci.lower' %in% names(res$coefficients)) {
            if ((!is.na(res$coefficients$z.ci.lower[cR])) & (!is.na(res$coefficients$z.ci.upper[cR]))) {
              tempfsqlw <- sprintf("%.1f", recursiveround(res$coefficients$z.ci.lower[cR], digits=1))
              if ((tempfsqlw == "-0.0") | (tempfsqlw == "0.0")) {
                tempfsqlw <- "0.0"
              }
              tempfsqup <- sprintf("%.1f", recursiveround(res$coefficients$z.ci.upper[cR], digits=1))
              if ((tempfsqup == "-0.0") | (tempfsqup == "0.0")) {
                tempfsqup <- "0.0"
              }
              temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                  floor(confidenceinterval*100), tempfsqlw, tempfsqup)
            }
          }
        }
      }
      
      temptext <- sprintf('%s, <span style="font-style: italic;">p</span>', temptext)
      temptext <- sprintf('%s %s %s', temptext, pullvalue$modifier, pullvalue$report)
      
      if (significanceconfidence) {
        if ('p.value.ci.lower' %in% names(res$coefficients)) {
          if ((!is.na(res$coefficients$p.value.ci.lower[cR])) & (!is.na(res$coefficients$p.value.ci.upper[cR]))) {
            tempfsqlw <- sprintf("%.3f", recursiveround(res$coefficients$p.value.ci.lower[cR], digits=3))
            tempfsqup <- sprintf("%.3f", recursiveround(res$coefficients$p.value.ci.upper[cR], digits=3))
            temptext <- sprintf('%s [%2.0f%% CI: %s to %s]', temptext,
                                floor(confidenceinterval*100), tempfsqlw, tempfsqup)
          }
        }
      }
      
      res$coefficients$text[cR] <- temptext
    }
    
  }
  
  modelsummary <- 'Regression analysis indicated that'
  if ('change' %in% names(res)) {
    modelsummary <- 'Hierarchical regression analysis indicated that'
  }
  if ('coefficients' %in% names(res)) {
    for (cR in 2:nrow(res$coefficients)) {
      tempstring <- res$coefficients$text[cR]
      tempstring <- stringr::str_split(tempstring, ", <span style=\"font-style: italic;\">t</span>")[[1]]
      modelsummary <- sprintf('%s %s (%s)', modelsummary, res$coefficients$Variable[cR], tempstring[1])
      if (cR < nrow(res$coefficients)) {
        if (cR < (nrow(res$coefficients)-1)) {
          modelsummary <- sprintf('%s, ', modelsummary)
        } else {
          modelsummary <- sprintf('%s and ', modelsummary)
        }
      }
    }
  } else {
    modelsummary <- sprintf('%s the model', modelsummary)
  }
  
  temppval <- Rmimic::fuzzyP(res$p.value)
  if ('change' %in% names(res)) {
    temppval <- Rmimic::fuzzyP(res$change$p.value)
  }
  outval <- paste(temppval$modifier,temppval$report,sep = " ")
  if (temppval$significance) {
    modelsummary <- sprintf("%s explained a statistically significant", modelsummary)
  } else {
    modelsummary <- sprintf("%s failed to explain a statistically significant", modelsummary)
  }
  tempstring <- res$text
  tempstring <- stringr::str_split(tempstring, ", <span style=\"font-style: italic;\">R")[[1]]
  modelsummary <- sprintf("%s (%s),", modelsummary, tempstring[1])
  
  if ('change' %in% names(res)) {
    modelsummary <- sprintf("%s change in variance in %s (", modelsummary, res$dv)
  } else {
    modelsummary <- sprintf("%s amount of variance in %s (", modelsummary, res$dv)
  }
  
  if ('change' %in% names(res)) {
    tempstring <- res$changetext
    tempstring <- stringr::str_split(tempstring, ", ")[[1]]
    tempstring <- tempstring[str_detect(tempstring, '<span style=\"font-style: italic;\">&Delta;R')]
  } else {
    tempstring <- res$text
    tempstring <- stringr::str_split(tempstring, ", ")[[1]]
    tempstring <- tempstring[str_detect(tempstring, '<span style=\"font-style: italic;\">R')]
  }
  for (cR in 1:length(tempstring)) {
    modelsummary <- sprintf("%s%s", modelsummary, tempstring[cR])
    if (cR == 1) {
      modelsummary <- sprintf('%s, ', modelsummary)
    }
  }
  modelsummary <- sprintf("%s)", modelsummary)
  
  if ('VIF' %in% names(res)) {
    if (res$VIF > 2.5) {
      modelsummary <- sprintf("%s \n[Warning: VIF is %.1f, suggesting a", modelsummary, recursiveround(res$VIF, digits=1))
      if (res$VIF > 2.5) {
        degree <- "slight"
      }
      if (res$VIF > 5) {
        degree <- "moderate"
      }
      if (res$VIF > 10) {
        degree <- "severe"
      }
      modelsummary <- sprintf("%s %s", modelsummary, degree)
      modelsummary <- sprintf("%s multicollinearity issue]", modelsummary)
    }
  }
  res$modelsummary <- sprintf("%s.", modelsummary)
  
  
  return(res)
}


recursiveround <- function(value, digits=2) {
  studywiseAlpharounding <- sort(seq(digits, (digits+3), by=1), decreasing=TRUE)
  for (cSWA in 1:length(studywiseAlpharounding)) {
    value <- round(value, digits=studywiseAlpharounding[cSWA])
  }
  return(value)
}

