#' descriptives
#'
#' @description Compute SPSS style descriptives.
#'
#' @param variables Variable name or list of variables to compute descriptives for.
#' @param groupvariable Variable name or list of variables to use for grouping.
#' @param data Data frame containing the variables of interest.
#' @param verbose Parameter to print all output to console. Default is TRUE.
#' @param verbosedescriptives Parameter to print descriptives output to console. Default is TRUE.
#' @param verbosefrequencies Parameter to print frequency output to console. Default is TRUE.
#'
#' @return
#' \item{datamatrix}{Data table of descriptive statistics.}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 3, 2020
#'
#' @importFrom stats shapiro.test complete.cases median sd
#' @importFrom Rmisc CI
#' 
#' @examples
#' 
#'     # Compute descriptive statistics
#'     desc <- descriptives(data=mtcars, groupvariable=FALSE, 
#'     verbose=FALSE, verbosedescriptives=FALSE,
#'     verbosefrequencies=FALSE)
#'
#'     # Compute descriptive statistics by group
#'     desc <- descriptives(data=mtcars, groupvariable=c("cyl"), 
#'     verbose=TRUE, verbosedescriptives=TRUE, 
#'     verbosefrequencies=FALSE)
#'
#' @export

descriptives <- function(data=FALSE, variables=FALSE, groupvariable=FALSE, verbose=TRUE, verbosedescriptives=TRUE, verbosefrequencies=TRUE) {

  if (verbose == FALSE) {
    verbosedescriptives <- FALSE
    verbosefrequencies <- FALSE
  }
  
  if (variables[1] == FALSE) {
    variables <- names(data)
    comparisondataframe <- as.data.frame(data)
  } else {
    compvariables <- variables
    if (groupvariable[1] != FALSE) {
      compvariables <- c(compvariables, groupvariable)
    }
    comparisondataframe <- data.frame(as.data.frame(data[,compvariables]))
    names(comparisondataframe) <- compvariables
  }
  
  # Determine grouping and nongrouping variables
  varlabels <- toupper(colnames(comparisondataframe))
  groupvariable <- toupper(groupvariable)
  groupingvariables <- list()
  nongroupingvariables <- list()
  for (cB in 1:(length(varlabels))) {
    if (length(base::intersect(varlabels[cB], groupvariable)) == 0) {
      nongroupingvariables <- c(nongroupingvariables, cB)
    } else {
      groupingvariables <- c(groupingvariables, cB)
    }
  }
  groupingvariables <- unlist(groupingvariables)
  groupvariable <- colnames(comparisondataframe)[groupingvariables]
  nongroupingvariables <- unlist(nongroupingvariables)
  nongroupvariables <- colnames(comparisondataframe)[nongroupingvariables]
  rm(varlabels)
  
  # Setup output dataframe
  if (length(groupingvariables) == 0) {
    datamatrix <- expand.grid(Variable = nongroupvariables)
    datamatrix$CollapsedName <- 'All'
  } else {
    # determine labels
    outcal <- 'datamatrix <- expand.grid('
    for (cB in 1:(length(groupingvariables))) {
      outcal <- sprintf('%s%s = levels(factor(as.character(comparisondataframe[,groupingvariables[%d]])))',outcal,groupvariable[cB],cB)
      if (cB < length(groupingvariables)) {
        outcal <- sprintf('%s, ', outcal)
      }
    }
    outcal <- sprintf('%s, Variable = nongroupvariables)', outcal)
    eval(parse(text=outcal))
    rm(outcal)
    datamatrix$CollapsedName <- NA
  }
  datamatrix$N <- NA
  datamatrix$Missing <- NA
  datamatrix$Mean <- NA
  datamatrix$Median <- NA
  datamatrix$SD <- NA
  datamatrix$SE <- NA
  datamatrix$Min <- NA
  datamatrix$Max <- NA
  datamatrix$Distribution <- NA
  datamatrix$Frequencies <- NA
  datamatrix$CI.Lower <- NA
  datamatrix$CI.Upper <- NA
  
  # subset data
  for (cB in 1:nrow(datamatrix)) {
    
    workingdf <- comparisondataframe
    
    # subset for grouping variables
    if (length(groupingvariables) != 0) {
      outcal <- ''
      for (cC in 1:length(groupingvariables)) {
        workingdf <- workingdf[which(workingdf[,groupvariable[cC]] == as.character(datamatrix[cB,cC])),]
        outcal <- sprintf('%s%s',outcal,datamatrix[cB,cC])
        if (cC < length(groupingvariables)) {
          outcal <- sprintf('%s:', outcal)
        }
      }
      rm(cC)
      datamatrix$CollapsedName[cB] <- outcal
    }
    
    # Extract only the data for the variable of interest
    subvect <- workingdf[,as.character(datamatrix$Variable[cB])]
    
    datamatrix$N[cB] <- length(subvect)
    datamatrix$Missing[cB] <- sum(is.na(subvect))
    
    if (is.numeric(subvect)) {
      subvect <- subvect[stats::complete.cases(subvect)]
      datamatrix$Mean[cB] <- mean(subvect)
      datamatrix$Median[cB] <- stats::median(subvect)
      datamatrix$SD[cB] <- stats::sd(subvect)
      datamatrix$SE[cB] <- stats::sd(subvect)/sqrt(length(subvect))
      if (length(subvect) > 1) {
        datamatrix$Min[cB] <- min(subvect)
        datamatrix$Max[cB] <- max(subvect)
      }
      
      tempci <- Rmisc::CI(subvect, ci = 0.95)
      datamatrix$CI.Lower[cB] <- tempci[3]
      datamatrix$CI.Upper[cB] <- tempci[1]
      
      testouttemp <- 1 # assume normal unless definitely not
      tryCatch(testouttemp <- stats::shapiro.test(subvect)$p.value, error=function(e){})
      if (testouttemp == 1) {
        datamatrix$Distribution[cB] <- NA
      } else {
        if (testouttemp > 0.05) {
          datamatrix$Distribution[cB] <- 'Normal'
        } else {
          datamatrix$Distribution[cB] <- 'Not Normal'
        }
      }
    }
    frequentcases <- sort(unique(unlist(as.character(subvect))))
    if ((length(frequentcases) > 0) & (length(frequentcases) < 10) & (length(subvect) > 0)) {
      freqoutcom <- ''
      for (cF in 1:length(frequentcases)) {
        freqoutcom <- sprintf('%s%s = %d', freqoutcom,frequentcases[cF], length(subvect[which(subvect == frequentcases[cF])]))
        if (cF < length(frequentcases)) {
          freqoutcom <- sprintf('%s, ', freqoutcom)
        }
      }
      datamatrix$Frequencies[cB] <- freqoutcom
      rm(cF, freqoutcom)
    }
    
    rm(workingdf, subvect)
  }
  rm(cB)
  
  if (verbose != FALSE) {
    
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
    
    # find out what needs to be printed
    grouplabels <- unique(as.character(datamatrix$CollapsedName))
    nongrouplabels <- unique(as.character(datamatrix$Variable))
    ngroups <- length(grouplabels)
    if (ngroups > 1) {
      if (sum(unlist(datamatrix$Missing), na.rm = TRUE) > 0) {
        vectnames <- c("Variable", "CollapsedName", "N", "Missing", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      } else {
        vectnames <- c("Variable", "CollapsedName", "N", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      }
      vectfreqnames <- c("Variable", "CollapsedName", "Class", "Count", "Percent")
    } else {
      if (sum(unlist(datamatrix$Missing), na.rm = TRUE) > 0) {
        vectnames <- c("Variable", "N", "Missing", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      } else {
        vectnames <- c("Variable", "N", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      }
      vectfreqnames <- c("Variable", "Class", "Count", "Percent")
    }
    sepgap <- matrix(floor(spansize/length(vectnames)), nrow=1, ncol=length(vectnames))
    colnames(sepgap) <- vectnames
    sepgap[1,1] <- sepgap[1,1] + 2
    
    
    if (verbosedescriptives != FALSE) {
      
      # Write header
      summarylabel <- "Descriptive Statistics"
      cat(sprintf("\n%s\n", summarylabel))
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      for (cC in 1:length(vectnames)) {
        if (vectnames[cC] == "CollapsedName") {
          cat(sprintf("%-*s",sepgap[1,cC],"Group"))
        } else {
          cat(sprintf("%-*s",sepgap[1,cC],vectnames[cC]))
        }
      }
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
      for (cB in 1:length(nongrouplabels)) {
        # Write Variable Label
        cat(sprintf("%s\n",nongrouplabels[cB]))
        for (cG in 1:length(grouplabels)) {
          submatrix <- datamatrix[which(datamatrix$Variable == nongrouplabels[cB]),] # subset for variable
          submatrix <- submatrix[which(submatrix$CollapsedName == grouplabels[cG]),] # subset for group
          cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:1]))), " "), collapse = "")))
          starval <- 2
          if (vectnames[2] == "CollapsedName") {
            # Indent properly
            cat(sprintf("%s\n",grouplabels[cG]))
            cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:2]))), " "), collapse = "")))
            starval <- 3
          }
          
          # Cycle Through Values
          for (cC in starval:length(vectnames)) {
            if ((is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
              # value is numeric and not missing
              pullvalue <- submatrix[1,sprintf("%s", vectnames[cC])]
              if (is.integer(pullvalue)) {
                pullvalue <- sprintf('%d', pullvalue)
              } else {
                pullvalue <- sprintf('%0.1f', round(pullvalue, digits=1))
              }
              cat(sprintf("%-*s",sepgap[1,cC],pullvalue))
            } else {
              if ((!is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
                cat(sprintf("%-*s",sepgap[1,cC],submatrix[1,sprintf("%s", vectnames[cC])]))
              } else {
                cat(sprintf("%-*s",sepgap[1,cC],"na"))
              }
            }
          }
          cat(sprintf("\n"))
        }
        if (cB < length(nongrouplabels)) {
          cat(sprintf("%s\n",paste(replicate(floor(spansize/3)-1, bigspancharacter), collapse = "")))
        }
      }
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
    }
    
    if (verbosefrequencies != FALSE) {
      # Check to see if we need to output frequencies
      boolfreq <- 0
      for (cC in 1:nrow(datamatrix)) {
        if (!is.na(datamatrix$Frequencies[cC])) {
          boolfreq <- 1
        }
      }
      if (boolfreq == 1) {
        freqdatamatrix <- datamatrix[which(!is.na(datamatrix$Frequencies)),]
        grouplabels <- unique(as.character(freqdatamatrix$CollapsedName))
        nongrouplabels <- unique(as.character(freqdatamatrix$Variable))
        ngroups <- length(grouplabels)
        
        sepgap <- matrix(sepgap[1,1:length(vectfreqnames)], nrow=1, ncol=length(vectfreqnames))
        colnames(sepgap) <- vectfreqnames
        
        # Write header
        spansize <- floor(spansize*.66)
        summarylabel <- "Frequency Table"
        cat(sprintf("\n%s\n", summarylabel))
        cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
        for (cC in 1:length(vectfreqnames)) {
          if (vectfreqnames[cC] == "CollapsedName") {
            cat(sprintf("%-*s",sepgap[1,3],"Group"))
          } else {
            cat(sprintf("%-*s",sepgap[1,3],vectfreqnames[cC]))
          }
        }
        cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
        
        for (cB in 1:length(nongrouplabels)) {
          # Write Variable Label
          cat(sprintf("%s\n",nongrouplabels[cB]))
          for (cG in 1:length(grouplabels)) {
            submatrix <- datamatrix[which(datamatrix$Variable == nongrouplabels[cB]),] # subset for variable
            submatrix <- submatrix[which(submatrix$CollapsedName == grouplabels[cG]),] # subset for group
            starval <- 1
            if (vectfreqnames[2] == "CollapsedName") {
              cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1]))), " "), collapse = "")))
              # Indent properly
              cat(sprintf("%s\n",grouplabels[cG]))
              starval <- 2
            }
            
            # Cycle Through Values
            if (!is.na(submatrix$Frequencies[1])) {
              freqvector <- unlist(strsplit(submatrix$Frequencies, ","))
              if (!is.na(freqvector[1])) {
                for (cC in 1:length(freqvector)) {
                  cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:starval]))), " "), collapse = "")))
                  
                  tempvect <- unlist(strsplit(freqvector[cC], "="))
                  cat(sprintf("%-*s\n",sepgap[1,3],gsub(" ", "", tempvect[1])))
                  cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:(starval+1)]))), " "), collapse = "")))
                  
                  pullvalue <- as.integer(gsub(" ", "", tempvect[2]))
                  pullvaluePerc <- sprintf('%0.1f', round((pullvalue/submatrix$N[1])*100, digits=1))
                  cat(sprintf("%-*s",sepgap[1,3],pullvalue))
                  cat(sprintf("%-*s",sepgap[1,3],pullvaluePerc))
                  if (cC < length(freqvector)) {
                    cat(sprintf("\n"))
                  }
                }
              } else {
                cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:starval]))), " "), collapse = "")))
                cat(sprintf("%s","No frequency output available"))
              }
              if (cG < length(grouplabels)) {
                cat(sprintf("\n"))
              }
            }
          }
          if (cB < length(nongrouplabels)) {
            cat(sprintf("\n%s\n",paste(replicate(floor(spansize/3)-1, bigspancharacter), collapse = "")))
          }
        }
        cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
        
      }
    }
  }
  
  return(datamatrix)
}


