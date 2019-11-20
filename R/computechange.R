#' computechange
#'
#' @description Function to compute the change or difference between conditions. 
#'
#' @param data Database containing data
#' @param method Parameter to select if the computation should be done as Change or PercentChange. Default is Change.
#' @param dependentvariable Dependent Variable label
#' @param subjectid Subject ID label.
#' @param within Within subjects factor labels to focus on.
#' @param relativeto Factor level to use as the baseline. Default is to use the lowest factor level.
#' @param control Column labels for factors to account for.
#'
#' @return
#' \item{data}{Database with an additional variable reflecitng the change in the dependent variable.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, November 20, 2019
#'
#' @examples
#'     mockdatabase <- data.frame("ID" = rep_len(1:20,60), "Time" = c(rep_len("Time1",20),rep_len("Time2",20),rep_len("Time3",20)), "X" = runif(60))
#'     mockdatabase <- computechange(mockdatabase, dependentvariable='X', subjectid='ID', within='Time')
#'
#' @export

computechange <- function(data, method=NULL, dependentvariable=NULL, subjectid=NULL, within=NULL, relativeto=NULL, control=NULL) {
  
  between <- NULL
  if ((!is.null(between)) | (!is.null(within))) {
    
    completedata <- data[,c(subjectid[1], dependentvariable[1], between[1], within[1], control)]
    completedata$completedata2datalink <- 1:nrow(completedata)
    
    # Assess what got fed into the function
    betweenvariableL <- length(between)
    withinvariableL <- length(within)
    if (!is.null(subjectid)) {
      completedata[,subjectid[1]] <- unlist(as.character(completedata[,subjectid[1]])) # force to characters
      indivparticipant <- sort(unique(as.character(completedata[,subjectid[1]])))
    } else {
      indivparticipant <- sort(unique(as.character(rownames(data))))
      if ((withinvariableL > 0) & (betweenvariableL == 0)) {
        # a within subjects parameter has been specified
        stop("Error: A within subjects factor has been requested, but no subjectid parameter was provided to ensure pairwise comparisons. Please provide a column of subject ids.")
      }
    }
    
    if (is.null(method)) {
      method <- 'Change'
    } else if (tolower(method) == tolower("PercentChange")) {
      if (withinvariableL == 0) {
        method <- 'PercentDifference' # between groups
      }
    }
    method <- tolower(method)
    
    # determine other factors to accountfor
    factorofinterest <- c(between[1], within[1])
    if (is.factor(completedata[,factorofinterest[1]])) {
      factorlevelsofinterest <- levels(completedata[,factorofinterest[1]])
    } else {
      factorlevelsofinterest <- sort(unique(unlist(as.character((completedata[,factorofinterest[1]])))))
    }
    if (is.null(relativeto)) {
      relativeto <- unlist(as.character(factorlevelsofinterest[1]))
    } else {
      if (!(relativeto %in% factorlevelsofinterest)) {
        stop(sprintf("Error: The factor level of interest %s does not appear to be present within %s.", relativeto, factorofinterest[1]))
      }
    }
    factorlevelsofinterest <- factorlevelsofinterest[which(factorlevelsofinterest != relativeto)]
    otherfactorsofinterest <- c(subjectid[1], control)
    
    # create output variable
    outvarname <- sprintf('%s.%sin%srelativeto%s', dependentvariable,method,factorofinterest[1],relativeto)
    tempcal <- sprintf('completedata$%s <- NA',outvarname)
    suppressWarnings(eval(parse(text=tempcal)))
    rm(tempcal)
    # Create comparison set
    basedata <- completedata[which(completedata[,factorofinterest[1]] == relativeto),] 
    basenames <- unlist(colnames(basedata))
    
    # Cycle through other levels
    for (cB in 1:length(factorlevelsofinterest)) {
      compdata <- completedata[which(completedata[,factorofinterest[1]] == factorlevelsofinterest[cB]),] 
      if (nrow(compdata) > 0) {
        # tweak names to facilitate merging
        compnames <- unlist(colnames(compdata))
        compnameslist <- setdiff(compnames, otherfactorsofinterest)
        for (cC in 1:length(compnameslist)) {
          compnames[which(compnames == compnameslist[cC])] <- sprintf('%sforcomparison', compnameslist[cC])
        }
        colnames(compdata) <- compnames
        tempdatbase <- merge(basedata,compdata, by=otherfactorsofinterest, all.x=TRUE, all.y=TRUE)
        for (cObs in 1:nrow(tempdatbase)) {
          pulloutputrow <- tempdatbase$completedata2datalinkforcomparison[cObs]
          pullbasedat <- tempdatbase[cObs,dependentvariable]
          pullcompdat <- tempdatbase[cObs,sprintf('%sforcomparison', dependentvariable)]
          if ((!is.na(pullbasedat)) & (!is.na(pullcompdat))) {
            tempoutdat <- NA
            if (tolower(method) == tolower("Change")) {
              tempoutdat <- pullcompdat - pullbasedat
            } else if (tolower(method) == tolower("PercentChange")) {
              tempoutdat <- ((pullcompdat - pullbasedat) / pullbasedat) * 100
            } else if (tolower(method) == tolower("PercentDifference")) {
              tempoutdat <- ((abs(pullcompdat - pullbasedat)) / ((pullcompdat + pullbasedat) / 2)) * 100
            }
            
            # send data to full database
            completedata[pulloutputrow,outvarname] <- tempoutdat
          }
        }
      }
    }
    
    tempcal <- sprintf('data$%s <- completedata$%s',outvarname, outvarname)
    suppressWarnings(eval(parse(text=tempcal)))
    rm(tempcal)
  }

  return(data)
}
