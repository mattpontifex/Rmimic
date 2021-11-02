#' pop_updateRmimic
#'
#' @description Popup window to update Rmimic
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, Nov 2, 2021
#' 
#' @importFrom pkgcond suppress_conditions
#' 
#' @export

pop_updateRmimic <- function() {
  pkgcond::suppress_conditions(invisible(tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})))
  pkgcond::suppress_conditions(invisible(tryCatch(devtools::install_github("mattpontifex/Rmimic", force=TRUE), error=function(e){})))
}