#' lmerPosthocCorrectionsubprocess
#'
#' @description Subprocess of lmerPosthoc by way of the subprocess lmerPosthocCorrection. Obtains data tables for all posthoc tests and their locations.
#'
#' @param res list containing data frames for stats, randomstats, and rsquared values
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, March 28, 2025
#' 
#' @export

lmerPosthocCorrectionsubprocess <- function(results) {
  
  colsofinterest <- c('p.value', 'significant', 'location', 'row', 'decomp', 'hold', 'contrast', 'idtag', 'textoutput')
  posthoclinkmatrix <- data.frame(matrix(NA, nrow=0, ncol=length(colsofinterest)))
  colnames(posthoclinkmatrix) <- colsofinterest
  
  anovacolsofinterest <- c('p.value', 'significant', 'location', 'row', 'factorsinvolved', 'idtag', 'textoutput')
  posthocanovamatrix <- data.frame(matrix(NA, nrow=0, ncol=length(anovacolsofinterest)))
  colnames(posthocanovamatrix) <- anovacolsofinterest
  
  if ('posthoc' %in% names(results)) {
    
    uniquenames <- names(results$posthoc)
    for (cUN in 1:length(uniquenames)) {
      
      tempelement <- results$posthoc[[uniquenames[cUN]]]
      
      if (is.data.frame(tempelement)) {
        
        # place data
        subposthoclinkmatrix <- data.frame(matrix(NA, nrow=nrow(tempelement), ncol=length(colsofinterest)))
        colnames(subposthoclinkmatrix) <- colsofinterest
        subposthoclinkmatrix$p.value <- tempelement$p.value
        subposthoclinkmatrix$significant <- tempelement$significant
        subposthoclinkmatrix$row <- seq(1,nrow(tempelement), by=1)
        subposthoclinkmatrix$location <- sprintf("%s", uniquenames[cUN])
        
        subposthoclinkmatrix$decomp <- tempelement$decomp
        subposthoclinkmatrix$hold <- tempelement$hold
        subposthoclinkmatrix$contrast <- tempelement$contrast
        subposthoclinkmatrix$idtag <- tempelement$idtag
        subposthoclinkmatrix$textoutput <- tempelement$textoutput
        
        # merge with main data
        posthoclinkmatrix <- rbind(posthoclinkmatrix, subposthoclinkmatrix)
        
      } else {
        # element is a structure
        
        if (length(tempelement) > 0) {
          
          # recursive call
          if ('posthoc' %in% names(tempelement)) {
            tempelement <- lmerPosthocCorrectionsubprocess(tempelement)
            tempelement$posthoclinkmatrix$location <- paste(uniquenames[cUN], tempelement$posthoclinkmatrix$location, sep='$posthoc$')
            
            # merge with main data
            posthoclinkmatrix <- rbind(posthoclinkmatrix, tempelement$posthoclinkmatrix)
          }
          
          # pull fixed effects anova stats
          tempelementanova <- tempelement$stats
          
          # place data
          subposthocanovamatrix <- data.frame(matrix(NA, nrow=nrow(tempelementanova), ncol=length(anovacolsofinterest)))
          colnames(subposthocanovamatrix) <- anovacolsofinterest
          subposthocanovamatrix$p.value <- tempelementanova$p.value
          subposthocanovamatrix$significant <- tempelementanova$significance
          subposthocanovamatrix$row <- seq(1,nrow(tempelementanova), by=1)
          subposthocanovamatrix$location <- sprintf("%s", uniquenames[cUN])
          
          subposthocanovamatrix$factorsinvolved <- tempelementanova$factorsinvolved
          subposthocanovamatrix$idtag <- tempelementanova$idtag
          subposthocanovamatrix$textoutput <- tempelementanova$textoutput
          
          # merge with main data
          posthocanovamatrix <- rbind(posthocanovamatrix, subposthocanovamatrix)
          
        }
      }
    }
    results$posthoclinkmatrix <- posthoclinkmatrix
    results$posthocanovamatrix <- posthocanovamatrix
  }
  return(results)
}