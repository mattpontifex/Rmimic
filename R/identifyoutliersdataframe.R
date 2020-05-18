#' identifyoutliersdataframe
#'
#' @description Function to identify outliers in a data frame based upon the interquartile range (as SPSS does for boxplots) and replace those values with NA. SPSS identifies mild outliers (1.5 x IQR) with circles and extreme outliers (3 x IQR) with asterisks.
#'
#' @param data Data frame containing the variables of interest.
#' @param variables Variable name or list of variables to impute.
#' @param iqrlimit The interquartile limit. Default is 3 consistent with SPSS' identification of extreme outliers.
#' @param direction Parameter to restrict identification of outliers to a particular tail. Default is two tailed. Parameters are: loweronly or upperonly.
#' @param minuniquesamples Minimum number of samples for outlier detection. Default is 10.
#' @param verbose Parameter to print the identified outliers.
#'
#' @return
#' \item{datta}{Data frame with outliers replaced with NA.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 18, 2020
#'
#' @export

identifyoutliersdataframe <- function(data, variables=FALSE, iqrlimit=3, direction="NA", minuniquesamples=10, verbose=TRUE) {
  
  if (variables[1] == FALSE) {
    variables <- names(data)
  }
  
  for (cVars in 1:length(variables)) {
    cColumns <- which(colnames(data) == variables[cVars])
    if (cColumns > 0) {
      datavector <- unlist(data[,cColumns])
      
      datavector <- Rmimic::identifyoutliers(datavector, iqrlimit=iqrlimit, direction=iqrlimit, minuniquesamples=minuniquesamples, verbose=verbose)
      
      data[,cColumns] <- datavector
    } # column name exists
  } # each variable
  
  return(data)
}
