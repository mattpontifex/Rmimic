#' decimalplaces
#'
#' @description Obtain the number of decimal places of precision in a vector
#'
#' @usage
#' decimalplaces(datavector)
#'
#' @param datavector Vector of data
#'
#' @return
#' \item{precisionvector}{Integer of maximum number of decimal places in a vector.}
#'
#' @examples
#'     precisioninteger <- decimalplaces(334.3410000000000000)
#'
#' @export

decimalplaces <- function(datavector) {
  
  datavector <- datavector[which(datavector != is.null(datavector))]
  precisioninteger <- 0
  for (cB in 1:(length(datavector))) {
  
    x = datavector[cB]
    tryCatch(temptprecisioninteger <- nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]]), error=function(e){temptprecisioninteger <- 0})
    tryCatch(if (temptprecisioninteger > precisioninteger) { precisioninteger <- temptprecisioninteger }, error=function(e){temptprecisioninteger <- 0})
    
  }
  return(precisioninteger)
}