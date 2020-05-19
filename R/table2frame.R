#' table2frame
#'
#' @description Helper function to blow up a table into a data frame.
#'
#' @param data Table with a column labeled freq.
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 19, 2020
#'
#' @export

table2frame <- function(data) {

  if (is.table(data)) {
    # blow it up
    temp <- as.data.frame(data)
    if (length(which(tolower(names(temp)) == (tolower('freq')))) > 0) {
      freqindx <- which(tolower(names(temp)) == (tolower('freq')))
      tempdata <- data.frame(matrix(NA,nrow=0,ncol=length(names(temp))))
      names(tempdata) <- names(temp)
      for (cR in 1:nrow(temp)) {
        if (as.integer(temp[cR,freqindx]) > 0) {
          # assuming that there is a count
          subtempdata <- data.frame(matrix(NA,nrow=as.integer(temp[cR,freqindx]),ncol=length(names(temp))))
          names(subtempdata) <- names(temp)
          for (cC in 1:ncol(temp)) {
            subtempdata[,cC] <- temp[cR,cC]
          }
          tempdata <- rbind(tempdata, subtempdata)
          rm(subtempdata)
        }
      }
      data <- tempdata[,which(tolower(names(tempdata)) != (tolower('freq')))]
      rm(tempdata)
    }
  }
  return(data)
}