#' lmerSimulate_conditionaldistribution
#'
#' @description Function to create a simulated dataset using simulate of the requested sample size.
#'
#' @param fit modelfit from the lmer function
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label. 
#' @param between text list specifying the between subjects variables
#' @param targetN numeric value specifying the number of subjects in the simulated dataset.
#' 
#' @return fit modelfit from the lmer function
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 13, 2025
#' 
#' @importFrom stats model.frame simulate
#'
#' @export

lmerSimulate_conditionaldistribution <- function(fit, dependentvariable=NULL, subjectid=NULL, between=NULL, targetN=NULL) {
  
  tempdbs <- stats::model.frame(fit)
  if (length(between) > 0) {
    tempdbs$newvarnameforBTW <- do.call(paste, c(tempdbs[between], sep="_by_"))
  } else {
    tempdbs$newvarnameforBTW <- 'one'
  }
  uniqueids <- unique(tempdbs[,subjectid])
  uniqueidsL <- length(uniqueids)
  uniqueBTW <- unique(tempdbs$newvarnameforBTW)
  uniqueBTWL <- rep_len(NA, length(uniqueBTW))
  for (cBTW in 1:length(uniqueBTW)) {
    currentcount <- 0
    for (cSID in 1:uniqueidsL) {
      checkindx <- which(tempdbs[,subjectid] == uniqueids[cSID] & tempdbs$newvarnameforBTW == uniqueBTW[cBTW])
      if (length(checkindx) > 0) {
        currentcount <- currentcount + 1
      }
    }
    uniqueBTWL[cBTW] <- currentcount
  }
  mainsmp <- NULL
  
  if (is.null(targetN)) {
    targetN <- uniqueBTWL
  } else {
    if (length(targetN) < length(uniqueBTWL)) {
      targetN <- rep_len(targetN[1], length(uniqueBTWL))
    }
  }
  # determine how many loops
  reps <- ceiling(max(targetN, na.rm=TRUE) / min(uniqueBTWL, na.rm=TRUE)) * 4 # increase so cases always need to be removed
  for (cReps in 1:reps) {
    smp <- stats::model.frame(fit)
    gg <- stats::simulate(fit,1)
    smp[,dependentvariable] <- gg[,1] # swap data in
    smp[,subjectid] <- paste(sprintf('IDCode_%d_', cReps), smp[,subjectid], sep='')  # generic ID
    
    if (is.null(mainsmp)) {
      mainsmp <- smp
    } else {
      mainsmp <- rbind(mainsmp, smp)
    }
  }
  
  # check totals
  if (length(between) > 0) {
    mainsmp$newvarnameforBTW <- do.call(paste, c(mainsmp[between], sep="_by_"))
  } else {
    mainsmp$newvarnameforBTW <- 'one'
  }
  continuebeyondloop <- FALSE
  while (!continuebeyondloop) {
    currentuniqueids <- unique(mainsmp[,subjectid])
    currentuniqueidsL <- length(currentuniqueids)
    currentuniqueBTWL <- rep_len(NA, length(uniqueBTW))
    for (cBTW in 1:length(uniqueBTW)) {
      currentcount <- 0
      for (cSID in 1:currentuniqueidsL) {
        checkindx <- which(mainsmp[,subjectid] == currentuniqueids[cSID] & mainsmp$newvarnameforBTW == uniqueBTW[cBTW])
        if (length(checkindx) > 0) {
          currentcount <- currentcount + 1
        }
      }
      currentuniqueBTWL[cBTW] <- currentcount
    }
    
    if (identical(targetN, currentuniqueBTWL)) {
      continuebeyondloop <- TRUE
    } else {
      for (cBTW in 1:length(uniqueBTW)) {
        if (targetN[cBTW] < currentuniqueBTWL[cBTW]) {
          # remove the needed number of cases
          casestoremove <- currentuniqueBTWL[cBTW] - targetN[cBTW]
          submainsmp <- mainsmp[which(mainsmp$newvarnameforBTW == uniqueBTW[cBTW]),]
          casetoremove <- submainsmp[sample(1:nrow(submainsmp), casestoremove), subjectid]
          mainsmp <- mainsmp[which(!(mainsmp[,subjectid] %in% casetoremove)),]
        }
      }
    }
  }
  mainsmp$newvarnameforBTW <- NULL # remove it
  
  # rerun model on new data
  newfit <- update(fit, data=mainsmp, evaluate = TRUE)
  
  return(newfit)
}