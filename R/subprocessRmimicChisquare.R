#' subprocessRmimicChisquare
#'
#' @description Subprocess to compute Chi-square analysis. For samples less than 1000, Fishers exact test statistic is used.
#'
#' @param data Data frame containing the variables of interest.
#' @param variables Variable name or list of variables to compute descriptives for.
#'
#' @return
#' \item{stats}{Summary of analysis.}
#' \item{tabularoutput}{Chi square table.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, March 6, 2020
#' 
#' @importFrom gmodels CrossTable
#' @importFrom stringr str_pad
#' @importFrom utils capture.output
#' @importFrom pkgcond suppress_conditions
#' 
#' @export

subprocessRmimicChisquare <- function(variables=FALSE, data=FALSE) {
  
  # debug variables
  #variables=c('treatment','improvement')
  #data=workingdata  
  #variables=c('carb','cyl')
  #data=mtcars
  #variables=c('MotherBMI','ChildBMI')
  #data=workingdata

  # Prepare output
  res <- list()
  
  if (variables[1] == FALSE) {
    variables <- names(data)
    cDF <- data
  } else {
    cDF <- data.frame(data[,variables])
    names(cDF) <- variables
  }
  cDF <- cDF[which(complete.cases(cDF)),]
  
  boolfish <- FALSE
  if (nrow(cDF) < 1000) {
    boolfish <- TRUE
  }
  results <- utils::capture.output(pkgcond::suppress_conditions(gmodels::CrossTable(cDF[,1], cDF[,2], fisher=boolfish, chisq=TRUE, expected=TRUE, sresid=TRUE, format="SPSS"), split = FALSE))
  
  stopindx <- which(pmatch(results, 'Statistics for All Table Factors') == 1)
  tabularoutput <- results[2:(stopindx-2)]
  
  # switch out labels for clarity
  swindx <- which(grepl("cDF[, 1]", tabularoutput, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE) == TRUE)
  tempvect <- tabularoutput[swindx]
  tempvectend <- gregexpr(pattern ='|',tempvect, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE)[[1]][1]
  tempvectend <- tempvectend - 2
  newlab <- stringr::str_pad(variables[1], width = tempvectend, side = "left")
  substr(tempvect, start = 1, stop = tempvectend+1) <- newlab[1:tempvectend]
  tabularoutput[swindx] <- tempvect
  
  tempvectend <- nchar(tempvect)
  swindx <- which(grepl("cDF[, 2] ", tabularoutput, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE) == TRUE)
  tempvect <- tabularoutput[swindx]
  tempvectstart <- gregexpr(pattern ='|',tempvect, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE)[[1]][1]
  tempvectstart <- tempvectstart + 1
  newlab <- stringr::str_pad(variables[2], width = tempvectend, side = "right")
  substr(tempvect, start = tempvectstart+1, stop = tempvectend) <- newlab
  tabularoutput[swindx] <- tempvect
  
  res$tabularoutput <- tabularoutput
  
  #snag chi-square and fishers stats
  dataframeout <- data.frame(matrix(NA,nrow=1,ncol=4))
  colnames(dataframeout) <- c('Chi-square', 'df', 'N', 'p.value')
  dataframeout[1,'N'] <- nrow(cDF)
  chiindx <- results[which(grepl("Pearson's Chi-squared test", results, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE) == TRUE) + 2]
  tempvect <- strsplit(chiindx, " ")[[1]]
  tempvect <- tempvect[which(tempvect != "")]
  tempvect <- tempvect[which(tempvect != "=")]
  dataframeout[1,'Chi-square'] <- as.numeric(tempvect[2])
  dataframeout[1,'df'] <- as.numeric(tempvect[4])
  dataframeout[1,'p.value'] <- as.numeric(tempvect[6])
  #http://www.biostathandbook.com/small.html
  if (nrow(cDF) < 1000) {
    fishindx <- FALSE
    tempcheck <- which(grepl("Alternative hypothesis: true odds ratio is not equal to 1" , results, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE) == TRUE)
    if (length(tempcheck) > 0) {
      fishindx <- results[tempcheck + 1]
    } else {
      tempcheck <- which(grepl("Alternative hypothesis: two.sided" , results, fixed=TRUE, ignore.case = FALSE, useBytes = TRUE) == TRUE)
      if (length(tempcheck) > 0) {
        fishindx <- results[tempcheck + 1]
      }
    }
    
    # switch out p values
    if (fishindx[1] != FALSE) {
      tempvect <- strsplit(fishindx, " ")[[1]]
      dataframeout[1,'p.value'] <- as.numeric(tempvect[length(tempvect)])
    }
  }
  temptext <- "Chi-square ("
  temptext <- sprintf('%s%d', temptext, dataframeout[1,'df'])
  temptext <- sprintf('%s, %d)', temptext, dataframeout[1,'N'])
  temptext <- sprintf('%s = %.1f', temptext, round(dataframeout[1,'Chi-square'], digits=1))
  outPvalue <- Rmimic::fuzzyP(as.numeric(dataframeout[1,'p.value']))
  temptext <- sprintf('%s, p %s %s', temptext, outPvalue$modifier, outPvalue$report)
  dataframeout$textoutput[1] <- temptext
  res$stats <- dataframeout
  
  return(res)
}


