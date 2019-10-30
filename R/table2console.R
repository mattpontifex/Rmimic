#' table2console
#'
#' @description Function that takes a data table and prints a nicely formatted table to the console
#'
#' @param dataframein Data frame to be printed
#' @param sepgap Data frame of the same size and with the same headers as the dataframein specifying the number of characters for each column.
#' @param spansize Number of characters on a line
#' @param headers Parameter to control if the variable labels should be printed. Default is TRUE.
#' @param alternate Parameter to control if the first column should have its own row. Default is FALSE.
#' @param seperators Parameter to control if the seperator lines should be shown. Default is FALSE.
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 29, 2019
#' 
#' @export

table2console <- function(dataframein, sepgap=NULL, spansize=95, headers=TRUE, alternate=FALSE, seperators=FALSE) {
  
  spancharacter <- "-"
  bigspancharacter <- " - "
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
    bigspancharacter <- " _ "
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
    bigspancharacter <- " - "
  }
  
  # determine spacing
  numberofcolumns <- length(colnames(dataframein))
  if (is.null(sepgap)) {
    sepgap <- data.frame(matrix(floor(spansize/numberofcolumns), nrow=1, ncol=numberofcolumns))
    # check spacing
    tempcheck <- matrix(NA, nrow=1, ncol=numberofcolumns)
    for (cC in 1:ncol(dataframein)) {
      tempcheck[,cC] <- max(nchar(dataframein[,cC]))
    }
    if (sum(tempcheck) < spansize) {
      paddsize <- floor((spansize-sum(tempcheck))/numberofcolumns)
      if (paddsize > 0) {
        for (cC in 1:ncol(dataframein)) {
          sepgap[1,cC] <- tempcheck[1,cC] + paddsize
        }
      } else {
        sepgap[1,] <- tempcheck[1,]
      }
    }
  }
  colnames(sepgap) <- colnames(dataframein)
  
  if (headers == TRUE) {
    # hard start line
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # write headers
    textvect <- Rmimic::cell2span(colnames(dataframein), sepgap=sepgap, spansize=spansize, trimwarn=FALSE)
    cat(sprintf('%s\n', textvect))
    
    # hard end line
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  # loop through each cell and pull data
  for (cR in 1:nrow(dataframein)) {
    if (alternate == FALSE) {
      textvect <- Rmimic::cell2span(dataframein[cR,], sepgap=sepgap, spansize=spansize, trimwarn=FALSE)
      cat(sprintf('%s\n', textvect))
    } else {
      # write label
      texttowrite <- unlist(as.character(dataframein[cR,1]))
      cat(sprintf('%s\n', texttowrite))
      # write data
      texttowrite <- unlist(as.character(dataframein[cR,2:ncol(dataframein)]))
      texttowrite <- c("empty", texttowrite)
      textvect <- Rmimic::cell2span(texttowrite, sepgap=sepgap, spansize=spansize, trimwarn=FALSE)
      cat(sprintf('%s\n', textvect))
    }
    if ((seperators == TRUE) & (cR < nrow(dataframein))) {
      cat(sprintf("%s\n",paste(replicate(spansize/3, bigspancharacter), collapse = "")))
    }
  }
  
  # write footer line
  cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  
}