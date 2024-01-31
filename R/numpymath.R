#' numpymath
#'
#' @description Function to apply mathmatical operations, vector1 operator vector2.
#'
#' @param vector1 data vector 1
#' @param operator mathmatical operator
#' @param vector2 data vector 2
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, September 13, 2021
#' 
#' @examples
#' 
#'   outputresult <- numpymath(PlantGrowth$weight, '*', 2.2)
#'   
#' @export
#' 
#' 
numpymath <- function(vector1, operator, vector2) {
  if (length(vector1) >= length(vector2)) {
    outvect <- rep(NA, length(vector1))
  } else {
    outvect <- rep(NA, length(vector2))
  }
  for (cR in 1:length(outvect)) {
    if (length(vector1) > 1) {
      if (tibble::is_tibble(vector1[cR])) {
        activevalue = vector1[[cR]]
      } else {
        activevalue = vector1[cR]
      }
    } else {
      activevalue = vector1
    }
    if (length(vector2) > 1) {
      if (tibble::is_tibble(vector2[cR])) {
        activevalue2 = vector2[[cR]]
      } else {
        activevalue2 = vector2[cR]
      }
    } else {
      activevalue2 = vector2
    }
    if ((is.numeric(activevalue)) && (!is.na(activevalue))) {
      if ((is.numeric(activevalue2)) && (!is.na(activevalue2))) {
        if (operator == '+') {
          outvect[cR] <- activevalue + activevalue2
        } else if (operator == '-') {
          outvect[cR] <- activevalue - activevalue2
        } else if (operator == '*') {
          outvect[cR] <- activevalue * activevalue2
        } else if (operator == '/') {
          outvect[cR] <- activevalue / activevalue2
        }
      }
    }
  }
  return(outvect)
}
