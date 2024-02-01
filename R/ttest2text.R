#' ttest2text
#'
#' @description Output t-test results in APA style format.
#'
#' @param datalist t-test list output
#' @param verbose Parameter to print all output to console. Default is TRUE.
#' 
#' @return A list with the following elements:
#' \item{statistics}{Data table of statistics.}
#' \item{text}{Statistics formatted in text.}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 7, 2019
#' 
#' @importFrom cli style_italic
#'
#' @export

ttest2text <- function(datalist, verbose=TRUE) {
  
  dataframeout <- data.frame(matrix(NA,nrow=1,ncol=14))
  colnames(dataframeout) <- c('statistic', 'df', 'p.value', 'conf.int.lower', 'conf.int.upper', 'alternative', 'method', 'z.value', 'effectsize', 'effectsize.conf.int.lower', 'effectsize.conf.int.upper', 'correlation', 'correlation.p.value', 'stud.conf.int')
  
  # Load stat
  puller <- FALSE
  tryCatch(dataframeout[1,'statistic'] <- datalist$statistic[1], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'df'] <- datalist$parameter[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'p.value'] <- datalist$p.value[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'conf.int.lower'] <- datalist$conf.int[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'conf.int.upper'] <- datalist$conf.int[[2]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'z.value'] <- datalist$z.value[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'effectsize'] <- datalist$effectsize[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'effectsize.conf.int.lower'] <- datalist$effectsize.conf.int.lower[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'effectsize.conf.int.upper'] <- datalist$effectsize.conf.int.upper[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'correlation'] <- datalist$correlation[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'correlation.p.value'] <- datalist$correlation.p.value[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'alternative'] <- datalist$alternative[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'method'] <- datalist$method[[1]], error=function(e){puller <- TRUE})
  tryCatch(dataframeout[1,'stud.conf.int'] <- datalist$stud.conf.int[[1]], error=function(e){puller <- TRUE})
  rm(puller)
  
  res <- list()
  res$statistics <- dataframeout
  res$text <- ""
  
  if (!is.na(dataframeout[[1,'method']])) {
    
    ttesttype <- "independent"
    ttestdist <- "parametric"
  
    if ((dataframeout[[1,'method']] == " Two Sample t-test") | (dataframeout[[1,'method']] == "Welch Two Sample t-test")) {
      #independent samples t-test
      ttesttype <- "independent"
      ttestdist <- "parametric"
    } else if (dataframeout[[1,'method']] == "Paired t-test") {
      #paired samples parametric test
      ttesttype <- "paired"
      ttestdist <- "parametric"
    } else if (dataframeout[[1,'method']] == "Wilcoxon rank sum test") {
      #mann-whitney t-test (independent samples)
      ttesttype <- "independent"
      ttestdist <- "nonparametric"
    } else if (dataframeout[[1,'method']] == "Wilcoxon signed rank test") {
      #Wilcoxon signed rank test (paired samples)
      ttesttype <- "paired"
      ttestdist <- "nonparametric"
    } else {
      ttesttype <- FALSE
      ttestdist <- FALSE
    }
    
    if (ttesttype != FALSE) {
    
      spancharacter <- "-"
      operatingsystem <- Sys.info()['sysname']
      if (operatingsystem == "Windows") {
        spancharacter <- "_"
      } else if (operatingsystem == "Darwin") {
        spancharacter <- "-"
      }
      
      # Setup initial output
      if (!is.na(dataframeout[[1,'statistic']])) {
        dataframeout[[1,'statistic']] <- sprintf('%.1f', round(abs(as.double(dataframeout[[1,'statistic']])), digits = 1))
      }
      if (!is.na(dataframeout[[1,'df']])) {
        dataframeout[[1,'df']] <- sprintf('%.1f', round(as.double(dataframeout[[1,'df']]), digits = 1)) 
        pullvalue <- dataframeout[[1,'df']]
        if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
          if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
            pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
          }
        }
        dataframeout[[1,'df']] <- pullvalue
      }
      if (ttestdist == "parametric") {
        res$text <- sprintf('%s(%s) = %s', cli::style_italic('t'), dataframeout[[1,'df']], dataframeout[[1,'statistic']])
      } else {
        if (ttesttype == "independent") {
          res$text <- sprintf('%s U = %s', cli::style_italic('Mann–Whitney U'), dataframeout[[1,'statistic']])
        } else {
          res$text <- sprintf('%s = %s', cli::style_italic('Wilcoxon V'), dataframeout[[1,'statistic']])
        }
        dataframeout[[1,'z.value']] <- sprintf('%.1f', round(as.double(dataframeout[[1,'z.value']]), digits = 1))
        res$text <- sprintf('%s, %s = %s', res$text, cli::style_italic('Z'), dataframeout[[1,'z.value']])
      }
      
      # report P value
      outPvalue <- fuzzyP(as.double(dataframeout[[1,'p.value']]))
      res$text <- sprintf('%s, %s %s %s', res$text, cli::style_italic('p'), outPvalue$modifier, outPvalue$report)

      # effect size reporting
      if (!is.na(dataframeout[[1,'effectsize']])) {
        if (is.na(dataframeout[[1,'stud.conf.int']])) {
          dataframeout[[1,'stud.conf.int']] <- '95%' # assume 95% confidence interval unless told otherwise
        } else {
          dataframeout[[1,'stud.conf.int']] <- sprintf('%2.0f%%', (round(as.double(dataframeout[[1,'stud.conf.int']]), digits = 2) * 100))
        }
        dataframeout[[1,'effectsize']] <- sprintf('%.2f', round(as.double(dataframeout[[1,'effectsize']]), digits = 2))
        if (ttestdist == "parametric") {
          if (ttesttype == "independent") {
            #if (operatingsystem == "Windows") {
            #  temptext <- sprintf("%s, d\u209b = %s", res$text, dataframeout[[1,'effectsize']])
            #} else {
            #  temptext <- sprintf("%s, ds = %s", res$text, dataframeout[[1,'effectsize']])
            #}
            temptext <- sprintf("%s, %s = %s", res$text, cli::style_italic(sprintf('%s%s', 'd','\u209b')), dataframeout[[1,'effectsize']])
            Encoding(temptext) <-  "UTF-8"
            res$text <- temptext
          } else {
            #if (operatingsystem == "Windows") {
            #  temptext <- sprintf('%s, d\u1d63\u2098 = %s', res$text, dataframeout[[1,'effectsize']])
            #} else {
            #  temptext <- sprintf('%s, drm = %s', res$text, dataframeout[[1,'effectsize']])
            #}
            temptext <- sprintf('%s, %s = %s', res$text, cli::style_italic(sprintf('%s%s%s', 'd','\u1d63', '\u2098')), dataframeout[[1,'effectsize']])
            Encoding(temptext) <-  "UTF-8"
            res$text <- temptext
          }
          if (!is.na(dataframeout[[1,'effectsize.conf.int.lower']])) {
            dataframeout[[1,'effectsize.conf.int.lower']] <- sprintf('%.2f', round(as.double(dataframeout[[1,'effectsize.conf.int.lower']]), digits = 2))
            dataframeout[[1,'effectsize.conf.int.upper']] <- sprintf('%.2f', round(as.double(dataframeout[[1,'effectsize.conf.int.upper']]), digits = 2))
            res$text <- sprintf('%s [%s CI: %s to %s].', res$text, dataframeout[[1,'stud.conf.int']], dataframeout[[1,'effectsize.conf.int.lower']], dataframeout[[1,'effectsize.conf.int.upper']])
          } else {
            res$text <- sprintf('%s.', res$text)
          }
        } else {
          # cohens r
          res$text <- sprintf('%s, %s = %s.', res$text, cli::style_italic('r'), dataframeout[[1,'effectsize']])
        }
      }
    }
    
    if (verbose != FALSE) {
      cat(sprintf('%s\n', res$text))
    }
  }
  
  return(res)
}  