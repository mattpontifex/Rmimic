#' antitibbler
#'
#' @description Function to return tibble to data frame.
#'
#' @param datain tibble of data
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, September 13, 2021
#'
#' @importFrom tibble is_tibble
#' @importFrom stringr str_replace_all
#' 
#' @export
#' 
#' 
antitibbler <- function(datain) {
  # check variable names for issues
  for (cR in 1:length(colnames(datain))) {
    tmpname <- colnames(datain)[cR]
    if ((tmpname == '') | (tmpname == sprintf('...%d', cR))) {
      # empty variable name fill with random characters
      tmpname <- sprintf('%s%04d%s', 
              paste0(sample(LETTERS, 3, TRUE), collapse=''),
              sample(9999, 1, TRUE), 
              paste0(sample(LETTERS, 3, TRUE), collapse=''))
    }
    if (Rmimic::fixname(tmpname) != colnames(datain)[cR]) {
      colnames(datain)[cR] <- Rmimic::fixname(tmpname)
    }
  }
  
  if (tibble::is_tibble(datain[1,1]) == TRUE) {
    dataout <- NA
    tmpcall <- sprintf('dataout <- data.frame(')
    for (cR in 1:length(colnames(datain))) {
      tmpcall <- sprintf("%s'%s' = datain$%s", tmpcall, colnames(datain)[cR], colnames(datain)[cR])
      if (cR < length(colnames(datain))) {
        tmpcall <- sprintf('%s,',tmpcall)
      }
    }
    tmpcall <- sprintf('%s)',tmpcall)
    eval(parse(text=tmpcall))
    
  } else {
    # data is not in tibble form
    dataout <- datain
  }
  return(dataout)
}