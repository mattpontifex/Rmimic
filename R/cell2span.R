#' cell2span
#'
#' @description Function that takes a vector of data table cells and creates a character string of a particular length. The text empty will result in an open span in the output
#'
#' @param dataframein Data frame to be transformed
#' @param sepgap Data frame of the same size and with the same headers as the dataframein specifying the number of characters for each column.
#' @param spansize Number of characters on a line
#' @param trimwarn Parameter to control if user should be warned that text was trimmed to fit. Default is TRUE
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 29, 2019
#' 
#' @export

cell2span <- function(dataframein, sepgap=NULL, spansize=95, trimwarn=TRUE) {
  
  # determine spacing
  if (is.data.frame(dataframein)) {
    dataframein <- as.vector(dataframein)
  }
  numberofcolumns <- length(dataframein)
    
  if (is.null(sepgap)) {
    sepgap <- data.frame(matrix(floor(spansize/numberofcolumns), nrow=1, ncol=numberofcolumns))
  }
  
  # populate output
  textvect <- sprintf("%s",paste(replicate(spansize, " "), collapse = ""))
  
  # fill output from dataframe
  warnissued <- FALSE
  for (cC in 1:numberofcolumns) {
    pullval <- unlist(as.character(dataframein[cC]))
    pullvalL <- nchar(pullval)
    if (pullval != "empty") {
      if (trimwarn == TRUE) {
        if ((pullvalL > sepgap[1,cC]) & (warnissued == FALSE)) {
          cat(sprintf('Warning from cell2span(): The cell contains data of length %d, but the requested output cell size was %d. Data was trimmed to fit the requested output size.\n', pullvalL, sepgap[1,cC]))
          warnissued <- TRUE
        }
      }
      # trim value to fit
      if (pullvalL > sepgap[1,cC]) {
        pullval <- substr(pullval, 1, sepgap[1,cC])
      }
      # pad value to fit
      if (pullvalL < sepgap[1,cC]) {
        # add trailing spaces
        numspacestoadd <- sepgap[1,cC] - pullvalL
        pullval <- sprintf('%s%s', pullval, paste(replicate(numspacestoadd, " "), collapse = ""))
      }
      
      # replace empty text with cell text
      startsp <- cC - 1
      if (cC == 1) {
        startsp <- 1
      } else {
        startsp <- sum(sepgap[1,1:cC-1])
      }
      stopsp <- sum(sepgap[1,1:cC])
      substr(textvect, startsp, stopsp) <- pullval
    }
  }
  
  return(textvect)
}