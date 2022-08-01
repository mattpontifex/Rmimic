#' ttesttopval
#'
#' @description P value for a t statistic.
#'
#' @param tval T test value.
#' @param df Degrees of Freedom.
#' @param tail Tail, default is two tailed. Options for Right or Left.
#' 
#' @return
#' \item{stats}{Summary of analysis.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 7, 2022
#'
#' @importFrom stats pt
#' 
#' @export

ttesttopval <- function(tval=FALSE,df=FALSE,tail=NULL) {  
  
  pval <- 1.0
  if (is.null(tail)) {
    pval <- 2*stats::pt(q=tval, df=df, lower.tail=FALSE)
  } else {
    if (tail == 'Right') {
      pval <- stats::pt(q=tval, df=df, lower.tail=FALSE)
    }
    if (tail == 'Left') {
      pval <- stats::pt(q=tval, df=df, lower.tail=TRUE)
    }
  }
  
  return(pval)
}