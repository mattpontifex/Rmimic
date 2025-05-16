#' lmerEffects_simpledesc
#'
#' @description Subfunction of Rmimic::lmerEffects. Creates descriptives table.
#'
#' @param results list output from the lmerEffects function containing the collection of results
#' 
#' @return results list output from the lmerEffects function
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom stats shapiro.test density.default model.frame
#' @importFrom purrr list_c
#'
#' @export

lmerEffects_simpledesc <- function(results) {
  
  tableout <- NULL
  
  workingdbs <- tryCatch({
    workingdbs <- results$fixedterms
  }, error = function(e) {
    workingdbs <- NULL
  })
  if (!is.null(workingdbs)) {
    # create basic table
    tempdbs <- stats::model.frame(results$fit)
    tempdbs <- tempdbs[,c(results$dependentvariable, results$fixedterms)]
    tempdbs$TempVarName149624052 <- NA
    for (cR in 1:nrow(tempdbs)) {
      tempdbs$TempVarName149624052[cR] <- paste(purrr::list_c(as.list(tempdbs[cR,2:(ncol(tempdbs)-1)])), collapse=" ")
    }
    tempdbs <- tempdbs[,c(results$dependentvariable, 'TempVarName149624052')]
    
    tableout <- data.frame(matrix(NA, nrow=length(unique(tempdbs$TempVarName149624052)), ncol=1))
    colnames(tableout) <- 'Group'
    tableout$Group <- unique(tempdbs$TempVarName149624052)
    
    tableout$N <- NA
    tableout$Missing <- NA
    tableout$Mean <- NA
    tableout$Median <- NA
    tableout$SD <- NA
    tableout$SE <- NA
    tableout$DistributionData <- NA
    tableout$DistributionDecision <- NA
    
    for (cR in 1:nrow(tableout)) {
      tempvect <- tempdbs[which(tempdbs$TempVarName149624052 == tableout$Group[cR]), 1] # DV is in column 1
      tableout$N[cR] <- length(tempvect)
      tableout$Missing[cR] <- sum(is.na(tempvect))
      tempvect <- tempvect[which(!is.na(tempvect))]
      
      tempout <- tryCatch({ tempout <- mean(tempvect, na.rm=TRUE) }, error = function(e) { tempout <- NULL })
      if (!is.null(tempout)) { tableout$Mean[cR] <- tempout }
      tempout <- tryCatch({ tempout <- median(tempvect, na.rm=TRUE) }, error = function(e) { tempout <- NULL })
      if (!is.null(tempout)) { tableout$Median[cR] <- tempout }
      tempout <- tryCatch({ tempout <- sd(tempvect, na.rm=TRUE) }, error = function(e) { tempout <- NULL })
      if (!is.null(tempout)) { tableout$SD[cR] <- tempout }
      if ((!is.na(tableout$SD[cR])) & !is.na(tableout$N[cR])) {
        tempout <- tryCatch({ tempout <- tableout$SD[cR] / sqrt(tableout$N[cR]) }, error = function(e) { tempout <- NULL })
        if (!is.null(tempout)) { tableout$SE[cR] <-tempout }
      }
      tempout <- tryCatch({ 
        # just give raw data
        #tempout <- paste(sort(tempdbs[which(tempdbs$TempVarName149624052 == tableout$Group[cR])],1), collapse = ",")
        densitycheck <- stats::density.default(tempdbs[which(tempdbs$TempVarName149624052 == tableout$Group[cR]),1], n=64, na.rm=TRUE)
        tempout <- paste(densitycheck$y, collapse = ",")
      }, error = function(e) { tempout <- NULL })
      if (!is.null(tempout)) { tableout$DistributionData[cR] <- tempout }
      
      # assume normal unless definitely not
      tempout <- tryCatch({ tempout <- stats::shapiro.test(tempvect)$p.value }, error = function(e) { tempout <- NULL })
      if (!is.null(tempout)) {
        if (tempout > 0.05) {
          tempout <- 'Normal'
        } else {
          tempout <- 'Not Normal'
        }
        tableout$DistributionDecision[cR] <- tempout 
      }
    }
    # dirty but effective
    tableout$Mean <- round(round(round(round(tableout$Mean, digits=4), digits=3), digits=2), digits=1)
    tableout$Median <- round(round(round(round(tableout$Median, digits=4), digits=3), digits=2), digits=1)
    tableout$SD <- round(round(round(round(tableout$SD, digits=4), digits=3), digits=2), digits=1)
    tableout$SE <- round(round(round(round(tableout$SE, digits=4), digits=3), digits=2), digits=1)
    
    tableout$Mean <- sprintf('%.1f', tableout$Mean)
    tableout$Median <- sprintf('%.1f', tableout$Median)
    tableout$SD <- sprintf('%.1f', tableout$SD)
    tableout$SE <- sprintf('%.1f', tableout$SE)
  }
  return(tableout)
}