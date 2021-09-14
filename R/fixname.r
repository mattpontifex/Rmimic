#' fixname
#'
#' @description Function to check variable names for issues and fix them.
#'
#' @param datain names
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, September 13, 2021
#'
#' @importFrom stringr str_replace_all
#' 
#' @export
#' 
#' 
fixname <- function(tmpname) {
  # check name
  problemcharacters <- c('!','@','#','$','%','^','&','*','.','<','>', ' ')
  for (cP in 1:length(problemcharacters)) {
    if (grepl(problemcharacters[cP], tmpname, fixed = TRUE) == TRUE) {
      tmpname <- stringr::str_replace_all(tmpname, problemcharacters[cP], '')
    }
  }
  # check for question mark and other special characters
  for (cP in 1:nchar(tmpname)) {
    if (substr(tmpname, start = cP, stop = cP) == '?') {
      tmpname <- stringr::str_replace_all(tmpname, '\\?', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == '(') {
      tmpname <- stringr::str_replace_all(tmpname, '\\(', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == ')') {
      tmpname <- stringr::str_replace_all(tmpname, '\\)', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == '[') {
      tmpname <- stringr::str_replace_all(tmpname, '\\[', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == ']') {
      tmpname <- stringr::str_replace_all(tmpname, '\\]', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == '{') {
      tmpname <- stringr::str_replace_all(tmpname, '\\{', '')
    }
    if (substr(tmpname, start = cP, stop = cP) == '}') {
      tmpname <- stringr::str_replace_all(tmpname, '\\}', '')
    }
  }
  return(tmpname)
}