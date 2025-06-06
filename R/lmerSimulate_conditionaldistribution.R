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
      smp <- data.frame(matrix(NA, nrow=0, ncol=ncol(tempdbs)))
      colnames(smp) <- colnames(tempdbs)
      for (cBTW in 1:length(uniqueBTW)) {
        btwsubtempdbs <- tempdbs[which(tempdbs$newvarnameforBTW == uniqueBTW[cBTW]),]
        btwuniqueids <- unique(btwsubtempdbs[,subjectid])
        btwuniqueidsL <- length(btwuniqueids)
        subuniqeids <- floor(btwuniqueidsL * subsample)
        selectsamples <- rep_len(TRUE, btwuniqueidsL)
        selectsamples[sample(1:btwuniqueidsL, btwuniqueidsL-subuniqeids, replace=FALSE)] <- FALSE
        includedids <- btwuniqueids[which(selectsamples)]
        btwsmp <- btwsubtempdbs[which(btwsubtempdbs[,subjectid] %in% includedids),]
        smp <- rbind(smp, btwsmp)
      }
        
      # Some piece of information is not getting updated - resulting in a downstream bug when we do posthoc analyses
      fit <- update(fit, data=smp, evaluate = TRUE)
      
      #fixedformula <- Reduce(paste, deparse(stats::formula(fit, fixed.only = TRUE)))
      #fixedformula <- stringr::str_split(stringr::str_remove_all(fixedformula, ' '), '~')[[1]]
      #randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
      #randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
      #textcall <- sprintf('fit <- pkgcond::suppress_conditions(lmerTest::lmer(formula = %s ~ %s + %s, data = smp))', dependentvariable[1], fixedformula[2], randomformula[2])
      #eval(parse(text=textcall))
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
  
  # remove excess participants - there should always be excess based upon the approach
  smp <- data.frame(matrix(NA, nrow=0, ncol=ncol(mainsmp)))
  colnames(smp) <- colnames(mainsmp)
  for (cBTW in 1:length(uniqueBTW)) {
    # remove the needed number of cases
    btwsubtempdbs <- mainsmp[which(mainsmp$newvarnameforBTW == uniqueBTW[cBTW]),]
    btwuniqueids <- unique(btwsubtempdbs[,subjectid])
    btwuniqueidsL <- length(btwuniqueids)
    selectsamples <- rep_len(TRUE, btwuniqueidsL)
    selectsamples[sample(1:btwuniqueidsL, btwuniqueidsL-targetN[cBTW], replace=FALSE)] <- FALSE
    includedids <- btwuniqueids[which(selectsamples)]
    btwsmp <- btwsubtempdbs[which(btwsubtempdbs[,subjectid] %in% includedids),]
    smp <- rbind(smp, btwsmp)
  }
  smp$newvarnameforBTW <- NULL # remove it
  
  # rerun model on new data
  #newfit <- update(fit, data=smp, evaluate = TRUE)
  
  # i think the dataframe is carrying forward additional information that is causing a random bug to occur
  # maybe not - but doesnt hurt to keep
  outputdatalist <- as.list(smp)
  outputdataframe <- as.data.frame(outputdatalist)
  
  return(outputdataframe)
}