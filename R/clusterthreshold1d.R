#' clusterthreshold1d
#'
#' @description Calculate contiguous clusters of locations in a 1D array that are above or below some threshold and with some minimum size.
#'
#' @usage
#' clusterthreshold1d(datain, crit = 0.05,
#'     clustersize = 5, direction = 'LessThan')
#'
#' @param datain Data vector.
#' @param crit Decimal representation of critical value. Default 0.05.
#' @param clustersize The cluster size of values meeting critial value.
#' @param direction Directional parameter for application of the critical value. Default is LessThan.
#'
#' @return A numeric logical vector corresponding to those points that meet the criteria.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, June 24, 2019
#'
#' @examples
#'     # Returns a vector corresponding to those values below 0.05 in 5 sample clusters
#'     result <- clusterthreshold1d(stats::runif(100), crit = 0.05,
#'         clustersize = 5, direction = 'LessThan')
#'
#' @export

clusterthreshold1d <- function(datain, crit = 0.05, clustersize = 5, direction = 'LessThan') {
  
  # Error Checks
  if (is.null(clustersize)) {
    stop(cat(sprintf('Error: clusterthreshold1d() - No clustersize specified.\n')))
  }
  if (!((is.null(nrow(datain))) | (is.null(ncol(datain))))) {
    stop(cat(sprintf('Error: clusterthreshold1d() - Input must be 1 dimensional.\n')))
  }
  # Force data to have 1 row and multiple columns  
  if (is.null(ncol(datain))) {
    datain <- t(datain) 
  }
  if (clustersize > ncol(datain)) {
    stop(cat(sprintf('Error: clusterthreshold1d() - Cluster size exceeds data size.\n')))
  }
  # Setup output
  res = matrix(0, nrow=nrow(datain), ncol=ncol(datain))
  # Find values exceeding criterion
  if (toupper(direction) == toupper('GreaterThan')) {
    threslist <- which(datain >= crit)
  } else {
    threslist <- which(datain <= crit)
  }
  # For each identified point
  for (thresi in 1:length(threslist)) {
    # Check that values can be computed
    thresselectstart <- threslist[thresi] - floor(clustersize/2)
    thresselectstop <- threslist[thresi] + floor(clustersize/2)
    if ((thresselectstart >= 1) & (thresselectstop <= length(datain))) {
      # Select surrounding data
      thresselect <- datain[thresselectstart:thresselectstop]
      if (toupper(direction) == toupper('GreaterThan')) {
        thresselect <- (thresselect >= crit)
      } else {
        thresselect <- (thresselect <= crit)
      }
      # See if all surrounding data meets the criteria
      if (any(is.na(thresselect))) {
        res[threslist[thresi]] <- 0
      } else {
        if (all(thresselect)) {
          res[threslist[thresi]] <- 1
        } else {
          res[threslist[thresi]] <- 0
        }
      }
    } else {
      # Unable to evaluate due to cluster size so default to no
      res[threslist[thresi]] <- 0
    }
  }
  rm(threslist, thresi, thresselectstart, thresselectstop, thresselect, datain, crit, clustersize, direction)
  return(res)
}
