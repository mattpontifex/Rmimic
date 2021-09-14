#' pseudobootstrapsamples
#'
#' @description For Exploratory Use Only. This function returns a dataframe with additional samples included drawn randomly from the existing database.
#'
#' @param data Data frame containing the variables of interest.
#' @param newsamples Integer representing the number of new samples to draw.
#' @param subjectid Subject ID label.
#' @param method Parameter - Add vs Replace - to determine if new samples are added to the existing dataframe or if the returned dataframe should be fully resampled.
#' @param seed Seed for the random number generator to obtain reproducible results. NA will turn this parameter off.
#'
#' @return
#' \item{database}{Database of variables with missing data replaced for the specified variables.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, November 20, 2020
#'
#' @importFrom pkgcond suppress_conditions
#' 
#' @examples
#'
#'     # Increase sample by 10 
#'     data <- pseudobootstrapsamples(data=airquality, newsamples=10, subjectid=NULL, method='Add', seed=NA)
#'
#' @export

pseudobootstrapsamples <- function(data, newsamples=10, subjectid=NULL, method=NULL, seed=NA) {
  
  debug <- FALSE
  if (debug) {
    data <- Rmimic::gradwakefulness
    newsamples <- 10
    subjectid <- NULL
    method <- "Replace"
    seed <- NA
    
  }
  
  if (is.null(method)) {
    method <- "Add"
  } else {
    if (toupper(method) == toupper("Add")) {
      method <- "Add"
    } else {
      method <- "Replace"
    }
  }
  
  if (is.data.frame(data)) {
    if (nrow(data) > 1) {
      
      # set random seed
      if (is.na(seed)) {
        set.seed(NULL)
      } else {
        set.seed(seed)
      }
      
      # check samples
      if (is.null(subjectid)) {
        uniquesamples <- 1:nrow(data)
      } else {
        uniquesamples <- unique(as.character(unlist(data[,subjectid])))
      }
      if (toupper(method) == toupper("Replace")) {
        if (is.null(subjectid)) {
          newsamples <- newsamples + nrow(data)
        } else {
          newsamples <- newsamples + length(uniquesamples)
        }
      }
      
      # populate list of new samples
      samplevect <- sample.int(length(uniquesamples), newsamples, replace = TRUE)
      
      if (toupper(method) == toupper("Add")) {
        # add
        if (is.null(subjectid)) {
          data <- rbind(data,data[samplevect,]) # just add the new data
        } else {
          # subjectid exists
          for (cP in 1:length(samplevect)) {
            tempdata <- data[which(data[subjectid] == uniquesamples[samplevect[cP]]),]
            data <- rbind(data, tempdata)
          }
        }
      } else {
        # replace
        if (is.null(subjectid)) {
          data <- data[samplevect,]
        } else {
          # subjectid exists
          origdata <- data
          for (cP in 1:length(samplevect)) {
            tempdata <- origdata[which(origdata[subjectid] == uniquesamples[samplevect[cP]]),]
            if (cP == 1) {
              data <- tempdata
            } else {
              data <- rbind(data, tempdata)
            }
          }
        }
      }
      
      return(data)
    }
  }
}

