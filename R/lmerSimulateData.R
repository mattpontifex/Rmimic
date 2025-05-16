#' lmerSimulateData
#'
#' @description Function to create a simulated dataset from an lmer model. When method is conditionaldistribution, the function simulates data by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters using the simulate function. The data is simulated multiple times and then subjects are removed until the desired sample sizes for each between subjects variable is obtained. When method is covariance, the function simulates data from a multivariate distribution using the covariance matrix using the MASS::mvrnorm function (when parametric is true) or using the mnonr::unonr function (when parametric is false).
#'
#' @param fit modelfit from the lmer function
#' @param between text list specifying the between subjects variables
#' @param within text list specifying any variables that should be considered as within subjects. If null assumes all factors are between subjects.
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label. If left NULL assumes all factors are between subjects
#' @param subsample decimal parameter 0 to 1 indicating the percent of the sample to use. default is the full sample. Ignored when method is conditionaldistribution
#' @param inflation decimal parameter indicating the group level multiplier for the final sample. default is no inflation
#' @param parametric boolean parameter indicating if the simulation should create a parametrically simulated dataset using the MASS::mvrnorm function (true) or a nonnormal simulated dataset using the mnonr::unonr function (false). Ignored when method is conditionaldistribution
#' @param method text specifying either conditionaldistribution (default) to use the simulate function or covariance to use the MASS::mvrnorm function or mnonr::unonr functions.
#'
#' @return results list output from the lmerEffects function
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 10, 2025
#'
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     newdata <- Rmimic::lmerSimulateData(fit, between=c('Group'), within=c('Time'), dependentvariable = 'Alertness', subjectid = 'PartID', inflation=1.0, method='conditionaldistribution')
#'     newdata2 <- Rmimic::lmerSimulateData(fit, between=c('Group'), within=c('Time'), dependentvariable = 'Alertness', subjectid = 'PartID', inflation=1.0, method='covariance', subsample=0.8)                
#'                
#' @importFrom dplyr slice_sample
#' @importFrom MASS mvrnorm
#' @importFrom stringr str_split
#' @importFrom mnonr mardia unonr
#' @importFrom stats model.frame formula cov runif
#' @importFrom data.table as.data.table
#'
#' @export

lmerSimulateData <- function(fit, between=NULL, within=NULL, dependentvariable=NULL, subjectid=NULL, subsample=NULL, inflation=NULL, parametric=TRUE, method=NULL) {
  
  # Each repetition of the function would create a new simulated dataset
  # this approach should downweight any individual subject
  # simulation also allows for growing the sample
  
  if (!is.null(method)) {
    if (toupper(method) == toupper("covariance")) {
      method = "covariance"
    } else if (toupper(method) == toupper("conditionaldistribution")) {
      method = "conditionaldistribution"
    } else {
      method = "conditionaldistribution"
    }
  } else {
    method = "conditionaldistribution"
  }
  
  if (method == "covariance") {
    # the dataset would be simulated based upon the covariance matrix of a subsample of data
    # covariance matrix is unique for each between subject level
    # random effects not including the participant id are considered as between subjects levels
    # covariance matrix is common for all within subject levels reflecting the inherent relationship of these levels
    
    # no longer using this approach but keeping the code
    # Extract model components
    #fixef_vals <- lme4::fixef(fit) # Fixed effects
    #ranef_vals <- lme4::ranef(fit) # Random effects
    #sigma_eps <- stats::sigma(fit)                # Residual SD
    tempdbs <- stats::model.frame(fit)
    subject_ids <- unique(tempdbs[,subjectid])
    randomfactor <- colnames(tempdbs)  
    randomfactor <- randomfactor[!(randomfactor %in% dependentvariable)]
    if (!is.null(within)) {
      within <- within[within %in% colnames(tempdbs)]
      randomfactor <- randomfactor[!(randomfactor %in% within)]
    }
    if (!is.null(between)) {
      between <- between[between %in% colnames(tempdbs)]
      randomfactor <- randomfactor[!(randomfactor %in% between)]
    }
    
    # Define number of subjects and observations per subject
    n_subjects <- length(subject_ids)
    #simplesubjectstats <- psych::statsBy(tempdbs, subjectid)
    #n_obs_per_subject <- mean(simplesubjectstats$n, na.rm=TRUE)  
    #sd_n_obs_per_subject <- sd(simplesubjectstats$n, na.rm=TRUE) + 0.01
    
    # for each between subjects factor
    checknames <- colnames(tempdbs) # treat everything as between
    checknames <- checknames[which(checknames != dependentvariable)] # remove DV
    checknames <- checknames[which(checknames != subjectid)] # remove subjectID
    if (!is.null(within)) {
      tempdbs$newvarnameforWT <- do.call(paste, c(tempdbs[within], sep="_by_"))
      checknames <- checknames[which(!(checknames %in% within))] # remove within
    } else {
      tempdbs$newvarnameforWT <- 'one'
    }
    if (length(between) > 0) {
      tempdbs$newvarnameforBTW <- do.call(paste, c(tempdbs[between], sep="_by_"))
      checknames <- checknames[which(!(checknames %in% between))] # remove between
    } else {
      tempdbs$newvarnameforBTW <- 'one'
    }
    uniqueBTW <- unique(tempdbs$newvarnameforBTW) 
    numberofsamples <- floor(n_subjects / length(uniqueBTW))
    if (!is.null(inflation)) {
      if (!is.numeric(inflation)) {
        inflation <- 1.0
      }
    }
      
    mastertempdbs <- tempdbs
    
    # apply random effects
    randomformula <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
    randomformula <- stringr::str_split(stringr::str_remove_all(randomformula, ' '), '~')[[1]]
    randomformula <- randomformula[which(randomformula != dependentvariable)]
    randomformula <- stringr::str_split(randomformula, "\\+|\\*|\\/")[[1]]
    randomformula <- randomformula[which(randomformula != '')]
    
    
    # create output database
    mastersim_data <- data.frame(matrix(NA, nrow=0, ncol=ncol(tempdbs)))
    colnames(mastersim_data) <- colnames(tempdbs)
    
    for (cBTW in 1:length(uniqueBTW)) {
      checkdata <- tempdbs[which(tempdbs$newvarnameforBTW == uniqueBTW[cBTW]),]
      uniqueWT <- unique(checkdata$newvarnameforWT)
      uniqueIDs <- unique(checkdata[,subjectid])
      numberofsamples <- floor(length(uniqueIDs))
      if (!is.null(inflation)) {
        numberofsamples <- floor(numberofsamples * inflation)
      }
      
      # create table for correlation matrix
      outmatrix <- data.frame(matrix(NA, nrow=length(uniqueIDs), ncol=length(uniqueWT)))
      colnames(outmatrix) <- uniqueWT
      outvarmatrix <- data.frame(matrix(NA, nrow=length(uniqueIDs), ncol=length(uniqueWT)))
      colnames(outvarmatrix) <- uniqueWT
      samplecount <- c()
      for (cP in 1:length(uniqueIDs)) {
        subsamplecount <- c()
        for (cWT in 1:length(uniqueWT)) {
          checkidx <- which(checkdata$newvarnameforWT == uniqueWT[cWT] & checkdata[,subjectid] == uniqueIDs[cP])
          if (length(checkidx) > 0) {
            tempout <- checkdata[checkidx, dependentvariable]
            outmatrix[cP, cWT] <- mean(tempout, na.rm=TRUE)
            subsamplecount <- c(subsamplecount, length(tempout))
            outvarmatrix[cP, cWT] <- sd(tempout, na.rm=TRUE)
          }
        }
        samplecount <- c(samplecount, mean(subsamplecount, na.rm=TRUE))
      }
      
      # if there is only a single sample the sd would have failed
      if (unique(samplecount) == 1) {
        for (cWT in 1:length(uniqueWT)) {
          checkidx <- which(checkdata$newvarnameforWT == uniqueWT[cWT])
          if (length(checkidx) > 0) {
            tempout <- checkdata[checkidx, dependentvariable]
            outvarmatrix[, cWT] <- sd(tempout, na.rm=TRUE)
          }
        }
      }
      
      # transform to variance
      outvarmatrix <- outvarmatrix^2
      for (cWT in 1:length(uniqueWT)) {
        tempvect <- outvarmatrix[,cWT]
        tempvectlimits <- quantile(tempvect, c(.10, .90), na.rm=TRUE)
        outvarmatrix[which(tempvect < tempvectlimits[1]),cWT] <- NA
        outvarmatrix[which(tempvect > tempvectlimits[2]),cWT] <- NA
      }
      
      if (!is.null(subsample)) {
        if (nrow(outmatrix) > 10) {
          if ((0.0 <= subsample) & (subsample <= 1.0)) {
            subsample <- 0.8
          }
          subnumberofsamples <- floor(nrow(outmatrix) * subsample)
        } else {
          subnumberofsamples <- nrow(outmatrix) - 1
        }
       # knock out 
        selectsamples <- rep_len(TRUE, nrow(outmatrix))
        selectsamples[sample(1:nrow(outmatrix), nrow(outmatrix)-subnumberofsamples, replace=FALSE)] <- FALSE
        outmatrix <- outmatrix[selectsamples,]
        outvarmatrix <- outvarmatrix[selectsamples,]
      }
      if (length(uniqueWT) > 1) {
        tempsds <- colMeans(outvarmatrix, na.rm=TRUE)
      } else {
        tempsds <- mean(outvarmatrix, na.rm=TRUE)
      }
      
      # grab the means and prevent an issue
      if (length(uniqueWT) > 1) {
        tempmeans <- colMeans(outmatrix, na.rm=TRUE)
        # stupid simple hack for insufficient observations
        outmatrix <- rbind(rbind(outmatrix, outmatrix), outmatrix)
        tempsigma <- stats::cov(outmatrix, use="pairwise.complete.obs")
      } else {
        tempmeans <- mean(outmatrix, na.rm=TRUE)
        # stupid simple hack for insufficient observations
        outmatrix <- rbind(rbind(outmatrix, outmatrix), outmatrix)
        tempsigma <- matrix(c(tempsds,1),1,1)
      }
      
      multiplecasespersubject <- mean(samplecount, na.rm=TRUE, trim=0.25)
      for (cMCpS in 1:multiplecasespersubject) {
        if (length(uniqueWT) > 1) {
          # Simulate Multivariate Data
          if (parametric) {
            withindata <- MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma)
          } else {
            skewkurtcheck <- mnonr::mardia(outmatrix)
            skewkurtcheckuni <- data.table::as.data.table(skewkurtcheck$univariate)
            rownames(skewkurtcheckuni) <- rownames(skewkurtcheck$univariate)
            withindata <- suppressWarnings(mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=unlist(skewkurtcheckuni[1,]), kurtosis = unlist(skewkurtcheckuni[3,])))
            colnames(withindata) <- colnames(outmatrix)
          }
          # unpack data
          for (cWT in 1:length(uniqueWT)) {
            submastersim_data <- data.frame(matrix(NA, nrow=numberofsamples, ncol=ncol(tempdbs)))
            colnames(submastersim_data) <- colnames(tempdbs)
            submastersim_data[,dependentvariable] <- withindata[,cWT]
            submastersim_data$newvarnameforBTW <- uniqueBTW[cBTW]
            submastersim_data$newvarnameforWT <- uniqueWT[cWT]
            submastersim_data[,subjectid] <- sprintf('IDCode_%s_%d', uniqueBTW[cBTW], seq(1:numberofsamples)) # generic ID
            # merge database
            mastersim_data <- rbind(mastersim_data, submastersim_data)
          }
        } else {
          # Simulate Univariate Data
          if (parametric) {
            withindata <- MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma)
          } else {
            withindata <- mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=NULL, kurtosis = NULL)
            colnames(withindata) <- colnames(outmatrix)
          }
          submastersim_data <- data.frame(matrix(NA, nrow=numberofsamples, ncol=ncol(tempdbs)))
          colnames(submastersim_data) <- colnames(tempdbs)
          submastersim_data[,dependentvariable] <- withindata[,1]
          submastersim_data$newvarnameforBTW <- uniqueBTW[cBTW]
          submastersim_data$newvarnameforWT <- uniqueWT[1]
          submastersim_data[,subjectid] <- sprintf('IDCode_%s_%d', uniqueBTW[cBTW], seq(1:numberofsamples)) # generic ID
          # merge database
          mastersim_data <- rbind(mastersim_data, submastersim_data)
        }
      } # cases per subject
    } # between subjects
    
    # decode
    if (length(between) > 0) {
      tempvect <- stringr::str_split(mastersim_data$newvarnameforBTW, '_by_')
      for (cN in 1:length(between)) {
        mastersim_data[,between[cN]] <- sapply(tempvect, function(x) { x[[cN]] })
      }
    }
    
    if (!is.null(within)) {
      tempvect <- stringr::str_split(mastersim_data$newvarnameforWT, '_by_')
      for (cN in 1:length(within)) {
        mastersim_data[,within[cN]] <- sapply(tempvect, function(x) { x[[cN]] })
      }
    }
    
    # restore unadjusted original data - just in case
    tempdbs <- mastertempdbs
    
    # populate other variables
    remainingfactors <- randomfactor
    remainingfactors <- remainingfactors[which(remainingfactors != subjectid)] # remove subjectID
    if (length(remainingfactors) > 0) {
      for (cRF in 1:length(remainingfactors)) {
        # does this factor ever vary within a participant
        withinfactorRF <- FALSE
        for (cUNs in 1:length(subject_ids)) {
          checkdata <- tempdbs[which(tempdbs[,subjectid] == subject_ids[cUNs]),]
          if (length(unique(checkdata[,remainingfactors[cRF]])) > 1) {
            withinfactorRF <- TRUE
          }
        }
        
        if (!(withinfactorRF)) {
          # between subjects factor
          
          # determine allocation ratio
          tempcal <- sprintf("subworkingdatabase <- doBy::summaryBy(%s ~ %s + %s, FUN=c(mean), data=tempdbs, keep.names=TRUE)", dependentvariable[1], subjectid[1], remainingfactors[cRF])
          suppressWarnings(eval(parse(text=tempcal)))
          
          uniqueRF <- unique(subworkingdatabase[,remainingfactors[cRF]])
          allocratio <- rep_len(NA, length(uniqueRF))
          for (cBTW in 1:length(uniqueRF)) {
            checkdata <- subworkingdatabase[which(subworkingdatabase[,remainingfactors[cRF]] == uniqueRF[cBTW]),]
            allocratio[cBTW] <- nrow(checkdata) / nrow(subworkingdatabase)
          }
          
          # each participant can only be 1
          totaloutsubjects <- length(unique(mastersim_data[,subjectid]))
          allocratio <- round(allocratio * totaloutsubjects, digits=0)
          if (sum(allocratio) < totaloutsubjects) {
            allocratio[1] <- allocratio[1] + 1
          }
          if (sum(allocratio) > totaloutsubjects) {
            allocratio[1] <- allocratio[1] - 1
          }
          # generate vector
          allocdistr <- c()
          for (cBTW in 1:length(uniqueRF)) {
            suballocdistr <- rep_len(sprintf('%s', uniqueRF[cBTW]), allocratio[cBTW])
            allocdistr <- c(allocdistr, suballocdistr)
          }
          allocdistr <- sample(allocdistr) # randomize
          
          uniqueoutsubjects <- unique(mastersim_data[,subjectid])
          for (cUNs in 1:length(uniqueoutsubjects)) {
            mastersim_data[which(mastersim_data[,subjectid] == uniqueoutsubjects[cUNs]), remainingfactors[cRF]] <- allocdistr[cUNs]
          }
        } # is a between subjects factor
      } # each remaining factor
    }
    
    applyrandomeffects <- FALSE
    # when not adding random - consistently best result
    #covariance: 100.0000%
    #original fixed: 13.1107% vs simmed: 14.4532%
    #original random: 27.7306% vs simmed: 28.0639%
    
    # when adding only random variances
    #covariance: 97.0000%
    #original fixed: 13.1107% vs simmed: 11.4930%
    #original random: 27.7306% vs simmed: 19.9025%
    
    # when adding only random means
    #covariance: 99.0000%
    #original fixed: 13.1107% vs simmed: 11.7508%
    #original random: 27.7306% vs simmed: 40.2483%
    
    # when adding both random variance and means
    #covariance: 93.0000%
    #original fixed: 13.1107% vs simmed: 9.9768%
    #original random: 27.7306% vs simmed: 32.4383%
    if (applyrandomeffects) {
      # apply random effects in
      for (cRF in 1:length(randomformula)) {
        tempstr <- stringr::str_replace_all(stringr::str_split(randomformula[cRF], "[|]")[[1]], '[(]|[)]', "")
        tempstr <- tempstr[which(tempstr != "1")]  
        tempdbs$newvarnameforRF <- do.call(paste, c(tempdbs[tempstr], sep="_by_"))
        mastersim_data$newvarnameforRF <- do.call(paste, c(mastersim_data[tempstr], sep="_by_"))
        
        uniqueRF <- unique(tempdbs$newvarnameforRF) 
        uniqueRFmaster <- unique(mastersim_data$newvarnameforRF) 
        if (length(uniqueRF) < 4) {
          # reflects a discrete factor
          
          for (cBTW in 1:length(uniqueRF)) {
            checkdata <- tempdbs[which(tempdbs$newvarnameforRF == uniqueRF[cBTW]),]
            # compute unique additional standard deviation
            tempsds <- sd(checkdata[,dependentvariable], na.rm=TRUE) - sd(tempdbs[,dependentvariable], na.rm=TRUE)
            # transform to variance
            tempsds <- tempsds^2
            # compute intercept as difference from overall mean
            tempmeans <- mean(checkdata[,dependentvariable], na.rm=TRUE) - mean(tempdbs[,dependentvariable], na.rm=TRUE)
            # mimic covariance matrix
            tempsigma <- matrix(c(tempsds,1),1,1)
            #tempmeans <- tempmeans * 0
            
            checkdataout <- mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]),]
            numberofsamples <- nrow(checkdataout)
            
            # Simulate Univariate Data
            if (parametric) {
              withindata <- MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma)
              withindata <- unlist(as.list(withindata))
            } else {
              withindata <- mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=NULL, kurtosis = NULL)
              withindata <- unlist(as.list(withindata))
            }
            
            # add the random effect in - only handles random intercepts - not slopes
            tempvect <- mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]), dependentvariable]
            mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]),dependentvariable] <- tempvect + withindata
            
          } #
        } else {
          # likely reflects some distribution of values
          potentialmeans <- c()
          potentialvariances <- c()
          
          # find the distribution
          for (cBTW in 1:length(uniqueRF)) {
            checkdata <- tempdbs[which(tempdbs$newvarnameforRF == uniqueRF[cBTW]),]
            # compute unique additional standard deviation
            tempsds <- sd(checkdata[,dependentvariable], na.rm=TRUE) - sd(tempdbs[,dependentvariable], na.rm=TRUE)
            # transform to variance
            tempsds <- tempsds^2
            # compute intercept as difference from overall mean
            tempmeans <- mean(checkdata[,dependentvariable], na.rm=TRUE) - mean(tempdbs[,dependentvariable], na.rm=TRUE)
            
            potentialmeans <- c(potentialmeans, tempmeans)
            potentialvariances <- c(potentialvariances, tempsds)
          }
          # generate limits based upon the original data
          potentialmeanlimits <- quantile(potentialmeans, c(.10, .90), na.rm=TRUE)
          #potentialmeanlimits <- c(min(potentialmeans, na.rm=TRUE), max(potentialmeans, na.rm=TRUE))
          potentialvariancelimits <- quantile(potentialvariances, c(.10, .90), na.rm=TRUE)
          #potentialvariancelimits <- c(min(potentialvariances, na.rm=TRUE), max(potentialvariances, na.rm=TRUE))
          
          for (cBTW in 1:length(uniqueRFmaster)) {
            # randomly generate a value from the potentials - all values are equally likely
            tempmeans <- stats::runif(1, potentialmeanlimits[1], potentialmeanlimits[2]+0.0001)
            tempsds <- stats::runif(1, potentialvariancelimits[1], potentialvariancelimits[2]+0.0001)
            #tempmeans <- tempmeans * 0
            
            # mimic covariance matrix
            tempsigma <- matrix(tempsds[1], nrow=1, ncol=1)
            checkdataout <- mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]),]
            numberofsamples <- nrow(checkdataout)
            
            # Simulate Univariate Data
            if (parametric) {
              withindata <- MASS::mvrnorm(numberofsamples, tempmeans, Sigma = tempsigma)
              withindata <- unlist(as.list(withindata))
            } else {
              withindata <- suppressWarnings(mnonr::unonr(n=numberofsamples, mu=tempmeans, Sigma = tempsigma, skewness=NULL, kurtosis = NULL))
              withindata <- unlist(as.list(withindata))
            }
            
            # add the random effect in - only handles random intercepts - not slopes
            tempvect <- mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]), dependentvariable]
            mastersim_data[which(mastersim_data$newvarnameforRF == uniqueRFmaster[cBTW]), dependentvariable] <- tempvect + withindata
            
          } # each factor in master
        } # number of factor levels
      } # each random element
    }
    
     # remove decode
     mastersim_data$newvarnameforBTW <- NULL # remove it
     mastersim_data$newvarnameforWT <- NULL # remove it 
     mastersim_data$newvarnameforRF <- NULL # remove it
     
    # return same order
    tempdbs <- stats::model.frame(fit)
    mastersim_data <- mastersim_data[,colnames(tempdbs)]
  } else {
    #method = "conditionaldistribution"
    
    if (is.null(inflation)) {
      inflation <- 1.0
    } else {
      if (!is.numeric(inflation)) {
        inflation <- 1.0
      }
    }
    
    tempdbs <- stats::model.frame(fit)
    if (length(between) > 0) {
      tempdbs$newvarnameforBTW <- do.call(paste, c(tempdbs[between], sep="_by_"))
    } else {
      tempdbs$newvarnameforBTW <- 'one'
    }
    uniqueids <- unique(tempdbs[,subjectid])
    uniqueBTW <- unique(tempdbs$newvarnameforBTW)
    uniqueBTWL <- rep_len(NA, length(uniqueBTW))
    for (cBTW in 1:length(uniqueBTW)) {
      currentcount <- 0
      for (cSID in 1:length(uniqueids)) {
        checkindx <- which(tempdbs[,subjectid] == uniqueids[cSID] & tempdbs$newvarnameforBTW == uniqueBTW[cBTW])
        if (length(checkindx) > 0) {
          currentcount <- currentcount + 1
        }
      }
      uniqueBTWL[cBTW] <- currentcount
    }
    numberofsamples <- floor(inflation * uniqueBTWL)
    
    newfit <- Rmimic::lmerSimulate_conditionaldistribution(fit, dependentvariable=dependentvariable, subjectid=subjectid, between=between, targetN=numberofsamples)
    mastersim_data <- stats::model.frame(newfit)
    
  }
  return(mastersim_data)
}





  
