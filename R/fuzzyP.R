#' fuzzyP
#'
#' @description Function to properly round p values. Returns several values allowing \cr the user to choose what rounded value is appropriate.
#'
#' @param datain Significance value.
#' @param alpha Critical value for alpha.
#' @param html Boolean parameter for output in html format.
#'
#' @return A list with the following elements:
#' \item{interpret}{Rounded p value for use in interpretations.}
#' \item{report}{Rounded p value for use in summary reports.}
#' \item{exact}{Exact p value.}
#' \item{modifier}{Equals sign unless the p value is less than .001.}
#' \item{significance}{Boolean for if the value is less than or equal to the critical alpha.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 28, 2017
#'
#' @examples
#'     # Round a randomly generated number between 0 and 1.
#'     outPvalue <- fuzzyP(stats::runif(1))
#'
#' @export

fuzzyP <- function(datain, studywiseAlpha=0.05, html=FALSE) {

  equalityparameters <- c('<', '=', '>')
  if (html) {
    equalityparameters <- c('&lt;', '&#61;', '&gt;')
  }
  decimalprecision <- nchar(strsplit(sub('0+$', '', as.character(studywiseAlpha)), ".", fixed = TRUE)[[1]][[2]])
  studywiseAlpharounding <- sort(seq(decimalprecision, (decimalprecision+3), by=1), decreasing=TRUE)
  
  modifier <- equalityparameters[2]
  if (!is.na(datain)) {
    if (is.numeric(datain)) {
      
      # round data based upon precision of alpha parameter
      checkvalue <- datain
      for (cSWA in 1:length(studywiseAlpharounding)) {
        checkvalue <- round(checkvalue, digits=studywiseAlpharounding[cSWA])
        if (cSWA == (length(studywiseAlpharounding)-1)) {
          if (studywiseAlpharounding[cSWA] == 3) {
            report <- sprintf('%0.3f', round(checkvalue, digits=studywiseAlpharounding[cSWA]))
          } else if (studywiseAlpharounding[cSWA] == 4) {
            report <- sprintf('%0.4f', round(checkvalue, digits=studywiseAlpharounding[cSWA]))
          } else if (studywiseAlpharounding[cSWA] == 5) {
            report <- sprintf('%0.5f', round(checkvalue, digits=studywiseAlpharounding[cSWA]))
          } else if (studywiseAlpharounding[cSWA] == 6) {
            report <- sprintf('%0.6f', round(checkvalue, digits=studywiseAlpharounding[cSWA]))
          }
        }
      }
      interpret <- checkvalue # rounded to 2 decimal places (if alpha = 0.05)
      if ((interpret > 0) & (interpret < 1)) {
        currentprecision <- tryCatch({
          currentprecision <- nchar(strsplit(sub('0+$', '', as.character(checkvalue)), ".", fixed = TRUE)[[1]][[2]])
        }, error = function(e) {
          currentprecision <- 2
        })
      } else {
        currentprecision <- 1
      }
      if (currentprecision > 5) {
        currentprecision <- 5
      }
      if (currentprecision == 0) {
        currentprecision <- 1
      }
      if (datain > 0.99) {
        report <- "0.9"
        modifier <- equalityparameters[3]
      } else if (datain > (studywiseAlpha*10)) {
        # doesnt really matter as it is not even close
        if (currentprecision == 2) {
          report <- sprintf('%0.1f', round(interpret, digits=(currentprecision-1)))
        } else if (currentprecision == 3) {
          report <- sprintf('%0.2f', round(interpret, digits=(currentprecision-1)))
        } else if (currentprecision == 4) {
          report <- sprintf('%0.3f', round(interpret, digits=(currentprecision-1)))
        } else if (currentprecision == 5) {
          report <- sprintf('%0.4f', round(interpret, digits=(currentprecision-1)))
        }
      } else {
        # round the close number as it is used for posthoc
        tempvect <- strsplit(sub('0+$', '', as.character(studywiseAlpha)), ".", fixed = TRUE)[[1]]
        lowerbound <- as.numeric(sprintf('%s.%s51', tempvect[1], tempvect[2]))
        if (as.numeric(report) < lowerbound) {
          if (currentprecision == 2) {
            interpret <- floor(as.numeric(report) * 100) / 100
          } else if (currentprecision == 3) {
            interpret <- floor(as.numeric(report) * 1000) / 1000
          } else if (currentprecision == 4) {
            interpret <- floor(as.numeric(report) * 10000) / 10000
          } else if (currentprecision == 5) {
            interpret <- floor(as.numeric(report) * 100000) / 100000
          }
        }
        # round
        lowerbound <- c(0.001, '0.001')
        if (currentprecision == 3) {
          lowerbound <- c(0.0001, '0.0001')
        } else if (currentprecision == 4) {
          lowerbound <- c(0.00001, '0.00001')
        } else if (currentprecision == 5) {
          lowerbound <- c(0.000001, '0.000001')
        }
        if (datain < lowerbound[1]) {
          report <- lowerbound[2]
          modifier <- equalityparameters[1]
        }
      }
    } else {
      report <- NA
      interpret <- 1
    }
  } else {
    report <- NA
    interpret <- 1
  }
  significance <- FALSE
  if (interpret <= studywiseAlpha) {
    significance <- TRUE
  }
  for (cInt in 1:3) {
    if (substr(report, nchar(report), nchar(report)) == "0") {
      report <- substr(report, 1, nchar(report)-1)
    }
  }
  
  res <- list()
  res$interpret <- interpret
  res$report <- report
  res$exact <- datain
  res$modifier <- modifier
  res$significance <- significance
  
  return(res)
}
