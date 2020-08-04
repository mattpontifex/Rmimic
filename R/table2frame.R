#' table2frame
#'
#' @description Helper function to blow up a table into a data frame.
#'
#' @param data Table with a column labeled freq.
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 19, 2020
#'
#' @examples
#' 
#'   # create frequency frame
#'   data <- data.frame("Outcome"="Yes","Predictor"="No","Freq"=3)
#'   data <- rbind(data, data.frame("Outcome"="No","Predictor"="No","Freq"=15))
#'   data <- rbind(data, data.frame("Outcome"="Yes","Predictor"="Yes","Freq"=20))
#'   data <- rbind(data, data.frame("Outcome"="No","Predictor"="Yes","Freq"=17))
#'   data <- table2frame(data)
#'   
#' @export

table2frame <- function(data) {
  # revised 8-4-2020 - added ability to handle data frame frequency table.

  #data <- data.frame("Outcome"="Yes","Predictor"="No","Freq"=3)
  #data <- rbind(data, data.frame("Outcome"="No","Predictor"="No","Freq"=15))
  #data <- rbind(data, data.frame("Outcome"="Yes","Predictor"="Yes","Freq"=0))
  #data <- rbind(data, data.frame("Outcome"="No","Predictor"="Yes","Freq"=17))
  
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
  } else if (is.data.frame(data)) {
    if (length(which(tolower(names(data)) == (tolower('freq')))) > 0) {
      freqindx <- which(tolower(names(data)) == (tolower('freq')))
      tempdata <- data.frame(matrix(NA,nrow=0,ncol=length(names(data))))
      names(tempdata) <- names(data)
      for (cR in 1:nrow(data)) {
        if (as.integer(data[cR,freqindx]) > 0) {
          # assuming that there is a count
          subtempdata <- data.frame(matrix(NA,nrow=as.integer(data[cR,freqindx]),ncol=length(names(data))))
          names(subtempdata) <- names(data)
          for (cC in 1:ncol(data)) {
            subtempdata[,cC] <- data[cR,cC]
          }
          tempdata <- rbind(tempdata, subtempdata)
          rm(subtempdata)
          
        } # if count is greater than zero
      } # loop through row
      
      data <- tempdata[,which(tolower(names(tempdata)) != (tolower('freq')))]
      rm(tempdata)
    } # has a column for freq
  } # is data frame
    
  return(data)
}