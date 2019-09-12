#' identifyoutliers
#'
#' @description Function to identify outliers based upon the interquartile range (as SPSS does for boxplots) and replace those values with NA. SPSS identifies mild outliers (1.5 x IQR) with circles and extreme outliers (3 x IQR) with asterisks.
#'
#' @param datavector Vector of data
#' @param iqrlimit The interquartile limit. Default is 3 consistent with SPSS' identification of extreme outliers.
#' @param direction Parameter to restrict identification of outliers to a particular tail. Default is two tailed. Parameters are: loweronly or upperonly.
#' @param minuniquesamples Minimum number of samples for outlier detection. Default is 10.
#' @param verbose Parameter to print the identified outliers.
#'
#' @return
#' \item{datavector}{Vector with outliers replaced with NA.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, September 12, 2019
#'
#' @importFrom stats quantile
#'
#' @examples
#'     # Populate generic data frame
#'     testvect <- runif(100,1,10)
#'
#'     # Insert outlying data points
#'     testvect[6] <- -1000
#'     testvect[10] <- 1000
#'
#'     testvect <- identifyoutliers(testvect, iqrlimit = 3, verbose=TRUE)
#'
#' @export

identifyoutliers <- function(datavector, iqrlimit=3, direction="NA", minuniquesamples=10, verbose=TRUE) {
  
  #library(utils)
  #tryCatch(library(stats), error=function(e){utils::install.packages("stats"); library(stats)}) # quantile
  
  
  # check if input is atomic
  if (is.atomic(datavector)) {
  
    # number of samples
    if (length(datavector) >= minuniquesamples) {
      
      # Calculate descriptive statistics
      temp <- stats::quantile(datavector, na.rm = TRUE, names = FALSE)
      tempIQR <- temp[4]-temp[2]
      
      # Determine outlying values
      upperoutlier <- temp[4] + (iqrlimit*tempIQR)
      loweroutlier <- temp[2] - (iqrlimit*tempIQR)
      upperout <- c()
      lowerout <- c()
      
      # Identify upper bounds
      if (toupper(direction) != toupper("loweronly")) {
        upperout <- which(datavector > upperoutlier)
        datavector[which(datavector > upperoutlier)] <- NA
      }
      
      # Identify lower bounds
      if (toupper(direction) != toupper("upperonly")) {
        lowerout <- which(datavector < loweroutlier)
        datavector[which(datavector < loweroutlier)] <- NA
      }
      
      if (verbose == TRUE) {
        cat(sprintf("identifyoutliers(): %.2f x IQR:\n", iqrlimit[1]))
        if (length(upperout)>0) {
          for (cR in 1:length(upperout)) {
            cat(sprintf("Case %d was an Upper Bound Outlier.\n", upperout[cR]))
          }
        }
        if (length(lowerout)>0) {
          for (cR in 1:length(lowerout)) {
            cat(sprintf("Case %d was an Lower Bound Outlier.\n", lowerout[cR]))
          }
        }
      }
    }
  }
  
  return(datavector)
}
