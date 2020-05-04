#' correlation
#'
#' @description Compute SPSS style correlation or partial correlation, with optional parameters for the approach.
#'
#' @param data Data frame containing the variables of interest.
#' @param variables Variable name or list of variables to compute descriptives for.
#' @param partial Parameter to specify the variable name to use for partial correlations.
#' @param method Specifies the type of correlation. Options are pearson (default), spearman, or kendall.
#' @param listwise Boolean operator for if the correlations should be computed using list-wise deletion for missing values.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param verbose Parameter to print all output to console. Default is TRUE.
#' 
#' @return A list with the following matrix elements:
#' \item{Comparison}{Comparisons made.}
#' \item{n}{Number of subjects.}
#' \item{estimate}{Correlation between comparisons.}
#' \item{\samp{p.value}}{Correlation p value.}
#' \item{statistic}{The value of the test statistic.}
#' \item{parameter}{The degrees of freedom of the test statistic.}
#' \item{\samp{conf.int.lower}}{Minimum confidence interval of correlation.}
#' \item{\samp{conf.int.upper}}{Maximum confidence interval of correlation.}
#' \item{method}{A character string indicating how the association was measured.}
#' \item{null.value}{The value of the association measure under the null hypothesis.}
#' \item{alternative}{A character string describing the alternative hypothesis.}
#' \item{textoutput}{A character string of the test result.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 4, 2020
#'
#' @importFrom stats cor.test
#' @importFrom ppcor pcor.test
#' @importFrom psych r.con
#'
#' @examples
#'     # Compute correlation.
#'     working <- data.frame("X" = runif(100), "Y" = runif(100), "Z" = runif(100))
#'     sR <- correlation(variables = c('X','Y','Z'), data = working,
#'         method="pearson", listwise=TRUE)
#'
#' @export

correlation <- function(variables=FALSE, partial=FALSE, data=FALSE, method=FALSE, listwise=TRUE, studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE) {
 
  # debug variables
  #data <- data.frame("X" = runif(100), "Y" = runif(100), "Z" = runif(100))
  #variables <- c('X','Y','Z')
  #partial <- FALSE
  #variables <- c('X','Y')
  #partial <- 'Z'
  #method <- FALSE
  #method <- "spearman"
  #method <- "kendall"
  #listwise <- TRUE
  #studywiseAlpha <- 0.05
  #confidenceinterval <- 0.95
  #verbose <- TRUE

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

  if (method[1] == FALSE) {
    method <- "pearson"
  } else if (toupper(method[1]) == toupper("spearman")) {
    method <- "spearman"
  } else if (toupper(method[1]) == toupper("kendall")) {
    method <- "kendall"
  }
  
  # prep database
  if (variables[1] == FALSE) {
    variables <- names(data)
    comparisondataframe <- data
  } else {
    compvariables <- variables
    if (partial[1] != FALSE) {
      compvariables <- c(compvariables, partial)
    }
    comparisondataframe <- data.frame(data[,compvariables])
    names(comparisondataframe) <- compvariables
  }
  varlabelsL <- length(variables)
  
  if (listwise == TRUE) {
    # remove NA from database
    comparisondataframe <- comparisondataframe[complete.cases(comparisondataframe),]
  }
  
  # create output structures
  res <- list()
  res <- data.frame(matrix(data=NA,nrow=1,ncol=15))
  names(res) <- c('Comparison','n','estimate', 'p.value','statistic','parameter','conf.int.lower','conf.int.upper','method','gp','null.value','alternative','comparison1','comparison2','textoutput')
  dataframeoutL <- 1
  
  # Run analysis
  for (cB in 1:varlabelsL) {
    for (cC in 1:varlabelsL) {
      if ((cB != cC) & (cB/cC >= 1)) {
        # add labels
        res[dataframeoutL,1] <- sprintf("%s and %s", variables[cC], variables[cB])
        if (partial[1]!=FALSE) {
          res[dataframeoutL,1] <- sprintf("%s: accounting for %s",res[dataframeoutL,1], partial)
        }
        res[dataframeoutL,'comparison1'] <- variables[cC]
        res[dataframeoutL,'comparison2'] <- variables[cB]
        
        # setup comparison
        compvariables <- c(variables[cC], variables[cB])
        if (partial[1] != FALSE) {
          compvariables <- c(compvariables, partial)
        }
        tempframe <- data.frame(comparisondataframe[,compvariables])
        tempframe <- tempframe[complete.cases(tempframe),]
        names(tempframe) <- compvariables
        res[dataframeoutL,2] <- nrow(tempframe)
        comparison1 <- tempframe[[1]]
        comparison2 <- tempframe[[2]]
        
        if (partial[1]==FALSE) {
          sR <- stats::cor.test(comparison2, comparison1, exact=FALSE, alternative='two.sided', method = method[1], conf.level = confidenceinterval, use = "complete.obs")
          if (method == "pearson") {
            res[dataframeoutL,'conf.int.lower'] <- sR$conf.int[1]
            res[dataframeoutL,'conf.int.upper'] <- sR$conf.int[2]
            res[dataframeoutL,'parameter'] <- sR$parameter[[1]]
          } else if (method == "spearman") {
            # Ruscio, J. (2008). Confidence intervals for Spearman's rank correlation with ordinal data:
            # A simulation study comparing analytic and bootstrap methods. Journal of Modern Applied
            # Statistical Analysis, 7(2), 416-434.
            res[dataframeoutL,'conf.int.lower'] <- tanh(((1/2)*log((1+sR$estimate)/(1-sR$estimate))) - sqrt(1/(nrow(tempframe)-3)) * qnorm((1+confidenceinterval)/2))
            res[dataframeoutL,'conf.int.upper'] <- tanh(((1/2)*log((1+sR$estimate)/(1-sR$estimate))) + sqrt(1/(nrow(tempframe)-3)) * qnorm((1+confidenceinterval)/2))
          }
          res[dataframeoutL,'gp'] <- 0
          res[dataframeoutL,'null.value'] <- sR$null.value[[1]]
          res[dataframeoutL,'alternative'] <- sR$alternative[[1]]
          res[dataframeoutL,'method'] <- sR$method[[1]]
        } else {
          comparison3 <- tempframe[,partial]
          sR <- ppcor::pcor.test(comparison2, comparison1, comparison3, method = method[1])
          if (method == "pearson") {
            temp <- psych::r.con(sR$estimate,nrow(tempframe),p=confidenceinterval,twotailed=TRUE)
            res[dataframeoutL,'conf.int.lower'] <- temp[1]
            res[dataframeoutL,'conf.int.upper'] <- temp[2]
            res[dataframeoutL,'method'] <- "Pearson's product-moment correlation"
          } else if (method == "spearman") {
            # Ruscio, J. (2008). Confidence intervals for Spearman's rank correlation with ordinal data:
            # A simulation study comparing analytic and bootstrap methods. Journal of Modern Applied
            # Statistical Analysis, 7(2), 416-434.
            res[dataframeoutL,'conf.int.lower'] <- tanh(((1/2)*log((1+sR$estimate)/(1-sR$estimate))) - sqrt(1/(nrow(tempframe)-3)) * qnorm((1+confidenceinterval)/2))
            res[dataframeoutL,'conf.int.upper'] <- tanh(((1/2)*log((1+sR$estimate)/(1-sR$estimate))) + sqrt(1/(nrow(tempframe)-3)) * qnorm((1+confidenceinterval)/2))
            res[dataframeoutL,'method'] <- "Spearman's rank correlation rho"
          } else if (method == "kendall") {
            res[dataframeoutL,'method'] <- "Kendall's rank correlation tau"
          }
          res[dataframeoutL,'gp'] <- sR$gp[[1]]
          res[dataframeoutL,'null.value'] <- 0
        }
        
        res[dataframeoutL,'estimate'] <- sR$estimate[[1]]
        res[dataframeoutL,'p.value'] <- sR$p.value[[1]]
        res[dataframeoutL,'statistic'] <- sR$statistic[[1]]
        
        # prepare p value for text output
        temppval <- Rmimic::fuzzyP(as.numeric(res[dataframeoutL,'p.value']))
        outval <- paste(temppval$modifier,temppval$interpret,sep = " ")
        
        # prepare strength statement
        if (abs(res[dataframeoutL,'estimate']) < 0.09) {
          interp <- "very weak"
        } else if (abs(res[dataframeoutL,'estimate']) < 0.26) {
          interp <- "weak"
        } else if (abs(res[dataframeoutL,'estimate']) < 0.39) {
          interp <- "weak to moderate"
        } else if (abs(res[dataframeoutL,'estimate']) < 0.6) {
          interp <- "moderate"
        } else if (abs(res[dataframeoutL,'estimate']) < 0.76) {
          interp <- "moderate to strong"
        } else if (abs(res[dataframeoutL,'estimate']) < 0.9) {
          interp <- "strong"
        } else if (abs(res[dataframeoutL,'estimate']) >= 0.9) {
          interp <- "very strong"
        }
          
        directmodifier <- ""
        if (res[dataframeoutL,'estimate'] < 0) {
          directmodifier <- " inverse"
        }
        phraser <- ""
        if (temppval$interpret <= studywiseAlpha) {
          phraser <- ""
        } else {
          phraser <- " non-significant"
        }
        
        if (toupper(method) == toupper("spearman")) {
          dtype <- "rs"
        } else if (toupper(method) == toupper("kendall")) {
          dtype <- "rt"
        } else {
          dtype <- "r"
        }
        
        # populate text output
        textout <- ''
        if (partial[1]!=FALSE) {
          textout <- sprintf('%sAfter accounting for variance associated with %s, there',textout,partial)
        } else {
          textout <- sprintf('%sThere',textout)
        }
        textout <- sprintf('%s was a %s%s%s',textout, interp, directmodifier, phraser)
        textout <- sprintf('%s relationship between %s and %s',textout,variables[cC], variables[cB])
        
        textout <- sprintf('%s, %s = %.2f',textout, dtype, res[dataframeoutL,'estimate'])
        if ((!is.na(res[dataframeoutL,'conf.int.lower'])) & (!is.na(res[dataframeoutL,'conf.int.upper']))) {
          textout <- sprintf('%s [%2.0f%% CI: %.2f to %.2f]',textout, confidenceinterval*100,
                             res[dataframeoutL,'conf.int.lower'], res[dataframeoutL,'conf.int.upper'])
        }
        textout <- sprintf('%s, p %s.',textout, outval)
        res[dataframeoutL,'textoutput'] <- textout
        
        # increment output level
        dataframeoutL <- dataframeoutL + 1
      }
    }
  }
  
  # output to console
  if (verbose == TRUE) {
    
    if (partial[1]!=FALSE) {
      temptext <- "Partial Correlations"
    } else {
      temptext <- "Correlations"
    }
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- "Analysis were conducted using"
    if (partial[1]==FALSE) {
      outstring <- sprintf('%s bivariate', outstring)
    }
    
    if (toupper(method) == toupper("spearman")) {
      outstring <- sprintf('%s Spearman\'s rho correlation coefficients', outstring)
    } else if (toupper(method) == toupper("kendall")) {
      outstring <- sprintf('%s Kendall\'s tau correlation coefficients', outstring)
    } else {
      outstring <- sprintf('%s Pearson\'s product moment correlation coefficients', outstring)
    }
    
    if (partial[1]!=FALSE) {
      outstring <- sprintf('%s accounting for the variance associated with %s', outstring, partial)
    }
    
    outstring <- sprintf('%s using the', outstring)
    if (partial[1]==FALSE) {
      outstring <- sprintf('%s stats (R Core Team, %s)', outstring, strsplit(as.character(utils::packageDate("stats")),"-")[[1]][1])
    } else {
      outstring <- sprintf('%s ppcor (Kim, %s)', outstring, strsplit(as.character(utils::packageDate("ppcor")),"-")[[1]][1])
      outstring <- sprintf('%s, psych (Revelle, %s),', outstring, strsplit(as.character(utils::packageDate("psych")),"-")[[1]][1])
    }
    outstring <- sprintf('%s and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    startindex <- 1
    segmentlength <- 8
    chunks <- ceiling(varlabelsL/segmentlength)
    for (cH in 1:chunks) {
      
      endindex <- startindex+segmentlength-1
      if (endindex > varlabelsL) {
        endindex <- varlabelsL
      }
        
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
      cat(sprintf("%-25s","Variable"))
      for (cB in startindex:endindex) {
        cat(sprintf("\t%5d.",cB))
      }
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
      # loop through each variable
      for (cB in 1:varlabelsL) {
        # output Variable name
        cat(sprintf("%d. %-44.40s\n",cB, variables[cB]))
        
        # output Correlation values
        cat(sprintf("%24s","Correlation"))
        for (cC in startindex:endindex) {
          if ((cB == cC) | (cB/cC < 1)) {
            cat(sprintf("\t%5s",""))
          } else {
            workingindx <- which(res$comparison1 == variables[cC] & res$comparison2 == variables[cB])
            tempval <- sprintf("%.3f", res$estimate[workingindx])
            temparray <- unlist(strsplit(as.character(tempval),"[.]"))
            tempval <- temparray[[2]]
            if (temparray[[1]] < 0) {
              tempval <- paste("-",".",tempval, sep = "")
            } else {
              tempval <- paste(".",tempval, sep = "")
            }
            cat(sprintf("\t%5s",tempval))
          }
        }
        cat(sprintf("\n"))
        
        # output p values
        cat(sprintf("%24s","Sig."))
        for (cC in startindex:endindex) {
          if (cB == cC) {
            cat(sprintf("\t%5s","----"))
          } else if (cB/cC < 1) {
            cat(sprintf("\t%5s",""))
          } else {
            workingindx <- which(res$comparison1 == variables[cC] & res$comparison2 == variables[cB])
            temppval <- fuzzyP(as.numeric(res$p.value[workingindx]))
            tempval <- sprintf("%.3f", round(temppval$exact,digits=3))
            temparray <- unlist(strsplit(as.character(tempval),"[.]"))
            outval <- paste(".",temparray[[2]], sep = "")
            if (temppval$interpret <= studywiseAlpha) {
              outval <- paste(outval,"*",sep = "")
            }
            if (temppval$modifier == "<") {
              tempval <- sprintf("%s", temppval$report)
              temparray <- unlist(strsplit(as.character(tempval),"[.]"))
              outval <- paste(".",temparray[[2]], sep = "")
              newval <- paste("<", outval, sep = "")
              outval <- newval
            }
            cat(sprintf("\t%5s",outval))
          }
        }
        cat(sprintf("\n"))
        
        if (listwise[1] == FALSE) {
          # output n values
          cat(sprintf("%24s","N"))
          for (cC in startindex:endindex) {
            if ((cB == cC) | (cB/cC < 1)) {
              cat(sprintf("\t%5s",""))
            } else {
              #cat(sprintf("\t%5s",Fn[cB,cC]))
              cat(sprintf("\t%5s",'99'))
            }
          }
          cat(sprintf("\n"))
        }
        cat(sprintf("\n"))
      }
      
      # Write Footer
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      cat(sprintf("* indicates correlation is significant at the %.2f level.\n", studywiseAlpha))
      if (listwise[1] != FALSE) {
        cat(sprintf("N = %d.\n",  nrow(comparisondataframe)))
      }
      
      if (chunks > 1) {
        cat(sprintf("\n\n"))
        startindex <- endindex+1
      } else {
        cat(sprintf("\n"))
      }
    } # end chunk
  
    outstring <- "Interpretation with Test Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    
    # check to see what gets reported
    testinterps <- c()
    for (cR in 1:nrow(res)) {
      temppval <- fuzzyP(as.numeric(res$p.value[cR]))
      testinterps <- c(testinterps, temppval$interpret)
    }
    if (sum(testinterps <= studywiseAlpha) > 0) {
      outstring <- "Statistically Significant Relationships"
      Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
      rm(outstring)
      Rmimic::typewriter(sprintf(paste(replicate(10, spancharacter), collapse = "")), tabs=2, spaces=0, characters=floor(spansize*.9))
      
      subres <- res[which(testinterps <= studywiseAlpha),]
      for (cR in 1:nrow(subres)) {
        Rmimic::typewriter(subres$textoutput[cR], tabs=2, spaces=0, characters=floor(spansize*.9))
        cat(sprintf("\n"))
      }
    }
    
    if (sum(testinterps > studywiseAlpha) > 0) {
      outstring <- "Non-significant Relationships"
      Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
      rm(outstring)
      Rmimic::typewriter(sprintf(paste(replicate(10, spancharacter), collapse = "")), tabs=2, spaces=0, characters=floor(spansize*.9))
      
      subres <- res[which(testinterps > studywiseAlpha),]
      for (cR in 1:nrow(subres)) {
        Rmimic::typewriter(subres$textoutput[cR], tabs=2, spaces=0, characters=floor(spansize*.9))
        cat(sprintf("\n"))
      }
    }
    
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  return(res)
}