#' fuzzyP
#'
#' @description Function to properly round p values. Returns several values allowing \cr the user to choose what rounded value is appropriate.
#'
#' @param datain Significance value.
#'
#' @return A list with the following elements:
#' \item{interpret}{Rounded p value for use in interpretations.}
#' \item{report}{Rounded p value for use in summary reports.}
#' \item{exact}{Exact p value.}
#' \item{modifier}{Equals sign unless the p value is less than .001.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 28, 2017
#'
#' @examples
#'     # Round a randomly generated number between 0 and 1.
#'     outPvalue <- fuzzyP(stats::runif(1))
#'
#' @export

fuzzyP <- function(datain) {

 modifier <- "="
 if (!is.na(datain)) {
   report <- sprintf('%0.3f', round(datain, digits=3))
   interpret <- report
    if (datain < 0.001) {
      interpret <- 0.001
      report <- '0.001'
      modifier <- "<"
    } else if (datain < 0.01) {
      interpret <- round(datain, digits=3)
    } else if (datain == 0.01) {
      interpret <- datain
    } else if (datain == 0.025) {
      interpret <- datain
    } else if (datain < 0.05) {
      interpret <- round(datain, digits=3)
    } else if (round(datain, digits=3) < 0.0551) {
      interpret <- floor(datain * 100) / 100
    } else if (datain < 0.06) {
      interpret <- round(datain, digits=3)
    } else if (datain < 0.1) {
      interpret <- round(datain, digits=2)
    } else if (datain < 0.7) {
      interpret <- round(datain, digits=1)
      report <- sprintf('%0.2f', round(datain, digits=2))
    } else {
      interpret <- round(datain, digits=1)
      report <- sprintf('%0.2f', round(datain, digits=1))
    }
   if (interpret == 1) {
     interpret <- round(datain, digits=2)
   }
   
   # remove trailing zeros in hundreds place
   if (substr(report, nchar(report), nchar(report)) == "0") {
     report <- substr(report, 1, nchar(report)-1)
   }
   # remove trailing zeros in tens place
   if (substr(report, nchar(report), nchar(report)) == "0") {
     report <- substr(report, 1, nchar(report)-1)
   }
   
 } else {
   report <- NA
   interpret <- 1
 }

 res <- list()
 res$interpret <- interpret
 res$report <- report
 res$exact <- datain
 res$modifier <- modifier
 class(res)

 return(res)
}
