#' multipleimputation
#'
#' @description Use multiple imputation (using the mice package) to replace missing data points. For each missing observation, multiple imputation is performed n times to create a combined multivariable modeling estimate. The database is returned with the missing data replaced with the modeling estimate.
#'
#' @param database Database of variables.
#' @param method A string specifying the imputation method. Default is predictive mean matching but see mice help file.
#' @param imputations The number of imputations to perform. Default is 50.
#' @param maxiterations The maximum number of iterations to perform. Default is 10 x the number of iterations.
#' @param restrict Parameter to restrict all estimates to the range of observed values.
#' @param seed Seed for the random number generator to obtain reproducible results. NA will turn this parameter off.
#'
#' @return
#' \item{database}{Database of variables with missing data replaced for the specified variables.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, August 18, 2017
#'
#' @import mice
#' @importFrom pkgcond suppress_conditions
#' 
#' @examples
#'
#'     # Run multiple imputation to replace missing cases
#'     # in built in dataset airquality
#'     data <- multipleimputation(airquality, imputations=10)
#'
#' @export

multipleimputation <- function(database, method="pmm", imputations=50, maxiterations=0, restrict=TRUE, seed=500) {

  # Use the mice package
  #Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained equations in R. Journal of Statistical Software, 45(3).
  #if (!require(mice)) { install.packages("mice"); library(mice) }

  if (maxiterations < 1) {
    maxiterations <- imputations*10
  }

  # Run multiple imputation
  tempData <- pkgcond::suppress_conditions(mice::mice(database, m=imputations, maxit=maxiterations, method=method, printFlag=FALSE, seed=seed))

  # Populate list of variable names
  for (cColumns in 1:(length(colnames(database)))) {

    if (tempData$nmis[[cColumns]] > 0) {
      
      # Obtain the precision
      temp <- tempData$data[,cColumns]
      precisioninteger <- NA
      if (is.numeric(temp)) {
        temp <- temp[complete.cases(temp)]
        precisioninteger <- Rmimic::decimalplaces(temp)
      }
      availablelevels <- sort(unique(unlist(as.character(temp))))
      
      # Loop through each missing datapoint
      missingpoints <- which(unlist(as.character(tempData$where[,cColumns])) == TRUE)
      for (cMissing in 1:length(missingpoints)) {
        
        # pull imputed results
        imputeddata <- tempData$imp[[cColumns]]
        imputeddata <- imputeddata[cMissing,]
        if (ncol(imputeddata) > 0) {
        
          outval <- NA
          if (!is.na(precisioninteger)) {
            # dealing with numbers
            
            # Restrict values to the range of observed values
            if (restrict == TRUE) {
              outvarmin <- min(temp, na.rm = TRUE)
              outvarmax <- max(temp, na.rm = TRUE)
              imputeddata[1,which(imputeddata[1,] < outvarmin)] <- NA
              imputeddata[1,which(imputeddata[1,] > outvarmax)] <- NA
            }
            
            # Compute the mean estimate
            outval <- as.numeric(sprintf("%.*f", precisioninteger, round(mean(as.numeric(imputeddata), na.rm = TRUE), digits = precisioninteger)))
            
          } else {
            # dealing with factors
            outval <- imputeddata[1]
          }
          
          # Replace the missing observation
          database[missingpoints[cMissing], cColumns] <- outval
        }
      }
    }
  }
  return(database)
}

