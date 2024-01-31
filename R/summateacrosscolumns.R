#' summateacrosscolumns
#'
#' @description Function to summate across columns
#'
#' @param workbook data frame that contains data
#' @param columnnames column names to sum across
#' @param makenumeric boolean parameter to force values to numeric
#' @param count boolean parameter to return the number of non empty values
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, September 13, 2021
#' 
#' @export
#' 
#' 
summateacrosscolumns <- function(workbook, columnnames, makenumeric=TRUE, count=FALSE) {
  outvect <- rep(NA, nrow(workbook))
  for (cR in 1:nrow(workbook)) {
    for (cC in 1:length(columnnames)) {
      if (tibble::is_tibble(workbook[cR,columnnames[cC]])) {
        activevalue = workbook[[cR,columnnames[cC]]]
      } else {
        activevalue = workbook[cR,columnnames[cC]]
      }
      if (makenumeric == TRUE) {
        activevalue <- unlist(as.numeric(activevalue))
      }
      if ((is.numeric(activevalue)) && (!is.na(activevalue))) {
        if (count == TRUE) {
          if (activevalue > 0) {
            activevalue <- 1
          } else {
            activevalue <- 0
          }
        }
        if (is.na(outvect[cR])) {
          outvect[cR] <- activevalue
        } else {
          outvect[cR] <- outvect[cR] + activevalue
        }
      } 
    }
  }
  return(outvect)
}