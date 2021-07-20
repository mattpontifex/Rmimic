#' dataframenamecheck
#'
#' @description Function to make sure dataframe names do not create an issue.
#'
#' @param datain dataframe
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 20, 2021
#' 
#' @importFrom stringr str_replace_all
#'
#'
#' @export

dataframenamecheck <- function(datain) {  
  
  tempnames <- names(datain)
  for (cV in 1:length(tempnames)) {
    tempnames[cV] <- stringr::str_replace_all(tempnames[cV], ' ', '')
    #tempnames[cV] <- stringr::str_replace_all(tempnames[cV], '.', '')
    #tempnames[cV] <- stringr::str_replace_all(tempnames[cV], '-', '')
    #tempnames[cV] <- stringr::str_replace_all(tempnames[cV], '_', '')
  }
  names(datain) <- tempnames
  
  return(datain)
}