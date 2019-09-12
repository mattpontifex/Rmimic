#' determineallpossiblecombinations
#'
#' @description Function to determine all possible factor combinations.
#'
#' @usage
#' determineallpossiblecombinations(datain)
#'
#' @param datain List of factors.
#'
#' @return A list of all possible combinations.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, November 13, 2017
#'
#' @examples
#'     # Determine all combinations of a generic list
#'     outlist <- determineallpossiblecombinations(c('A', 'B', 'C'))
#'
#' @export

determineallpossiblecombinations <- function(datain) {
  
  fixedterms <- as.character(datain) # force to character
  fullfixedterms <- c(fixedterms)
  if (length(fixedterms) > 1) {
    # Loop through potential combinations
    for (cL in 1:(length(fixedterms)-1)) {
      fixedtermsinteractions <- c()
      # Loop through all fixed terms
      for (cT1 in 1:length(fixedterms)) {
        # Loop through all created terms
        for (cT2 in 1:length(fullfixedterms)) {
          tempvect <- paste(sort(unique(c(fixedterms[cT1],fullfixedterms[cT2]))),collapse=":") # push them together
          tempvect <- unique(unlist(strsplit(as.character(tempvect),":"))) # split them back up, remove elements that are not unique
          fixedtermsinteractions <- unique(c(fixedtermsinteractions,sprintf("%s",paste(sort(tempvect), collapse=":"))))
        }
      }
      fullfixedterms <- unique(c(fullfixedterms,fixedtermsinteractions)) # remove nonunique factors
    }
  }
  
  # Check arrangement of variables in each element 
  for (cT2 in 1:length(fullfixedterms)) {
   tempvect <- unlist(strsplit(as.character(fullfixedterms[cT2]),":")) # split up element
   if (length(tempvect) > 1) {
     newtempvect <- c()
     for (cT1 in 1:length(fixedterms)) { # loop through each fixed term and extract them in order
       tempindex <- which(tempvect == fixedterms[cT1])
       if (length(tempindex) > 0) {
         newtempvect <- c(newtempvect,tempvect[tempindex])
       }
     }
     fullfixedterms[cT2] <- paste(newtempvect, collapse=":")
   }
  }
  return(fullfixedterms)
}

