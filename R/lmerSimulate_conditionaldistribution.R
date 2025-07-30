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
#' @importFrom data.table data.table as.data.table
#'
#' @export

lmerSimulate_conditionaldistribution <- function(fit, dependentvariable=NULL, subjectid=NULL, between=NULL, targetN=NULL, subsample=NULL) {
  
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
  

  if (!is.null(subsample)) {
    if (subsample < 1.0) {
      
      # subsample participants within each between subjects group
      tempdbs <- data.table::as.data.table(tempdbs)
      smplist <- vector("list", length(uniqueBTW))
      for (cBTW in seq_along(uniqueBTW)) {
        btwsubtempdbs <- tempdbs[newvarnameforBTW == uniqueBTW[cBTW]]
        btwuniqueids <- unique(btwsubtempdbs[[subjectid]])
        btwuniqueidsL <- length(btwuniqueids)
        subuniqeids <- floor(btwuniqueidsL * subsample)
        includedids <- sample(btwuniqueids, subuniqeids)
        btwsmp <- btwsubtempdbs[btwsubtempdbs[[subjectid]] %in% includedids]
        smplist[[cBTW]] <- btwsmp
      }
      smp <- data.table::rbindlist(smplist, use.names = TRUE, fill = TRUE)
      
      fit <- update(fit, data=smp, evaluate = TRUE)
    }
  }
  
  # determine how many loops
  reps <- ceiling(max(targetN, na.rm=TRUE) / min(uniqueBTWL, na.rm=TRUE)) * 2 # increase so cases always need to be removed
  
  simlist <- vector("list", reps)
  for (cReps in seq_len(reps)) {
    smp <- data.table::as.data.table(stats::model.frame(fit))
    gg <- stats::simulate(fit, 1)
    smp[, (dependentvariable) := gg[[1]] ]
    smp[, (subjectid) := paste(sprintf('IDCode_%d_', cReps), get(subjectid), sep='')] # generic ID
    simlist[[cReps]] <- smp
  }
  mainsmp <- data.table::rbindlist(simlist)
  
  # check totals
  if (length(between) > 0) {
    mainsmp[, newvarnameforBTW := do.call(paste, c(.SD, sep = "_by_")), .SDcols = between]
  } else {
    mainsmp[, newvarnameforBTW := "one"]
  }
  
  # remove excess participants - there should always be excess based upon the approach
  mainsmp <- data.table::as.data.table(mainsmp)
  smplist <- vector("list", length(uniqueBTW))
  for (cBTW in seq_along(uniqueBTW)) {
    N <- targetN[cBTW]
    btwsubtempdbs <- mainsmp[newvarnameforBTW == uniqueBTW[cBTW]]
    btwuniqueids <- unique(btwsubtempdbs[[subjectid]])
    
    includedids <- sample(btwuniqueids, N)
    btwsmp <- btwsubtempdbs[get(subjectid) %in% includedids]
    smplist[[cBTW]] <- btwsmp
  }
  smp <- data.table::rbindlist(smplist, use.names = TRUE)
  smp[, newvarnameforBTW := NULL]  # Remove column
  
  # return same order
  tempdbs <- stats::model.frame(fit)
  data.table::setcolorder(smp, colnames(tempdbs))
  smp <- smp[, .SD, .SDcols = colnames(tempdbs)]
  smp <- as.data.frame(smp)
  
  return(smp)
}