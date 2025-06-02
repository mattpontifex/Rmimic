#' colorpalette
#'
#' @description Function to return a set of standard color pallets
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, June 2, 2025
#'
#' @importFrom grDevices colorRampPalette
#' 
#' @examples
#'     # Round a randomly generated number between 0 and 1.
#'     workingpallete <- Rmimic::colorpalette()
#'
#' @export

colorpalette <- function() {
  
  res <- list()
  res$mainpallet <- c('#748FAC','#78B86E', 'gray50', '#746B98') 
  res$subpallet <- c('#677F99', '#6CA663', 'gray50', '#6B628C') 
  res$control <- '#748FAC'
  res$active <- '#78B86E'
  res$controlsub <- '#677F99'
  res$activesub <- '#6CA663'
  
  res$alt <- c('#540d6e', '#ee4266', '#ffd23f', '#3bceac')
  res$altsub <- c('#746B98','#EF6F63','#d9b552', '#1E9B8A') 
  
  res$fitline <- "#5c5c5c"
  res$fitfill <- "#aba9a9"
  
  res$crushparula <- grDevices::colorRampPalette(c("#00004B", "#1C2C75",'#38598C',	'#2B798B', '#1E9B8A',	'#85D54A',	'#FDE725', '#F9FB0E'))
  res$crushparulaalpha <- grDevices::colorRampPalette(c("#00004B50", "#1C2C7560",'#38598C70',	'#2B798B80', '#1E9B8A90',	'#85D54A',	'#FDE725', '#F9FB0E'))
  res$greenstain <- grDevices::colorRampPalette(c('#007f5f', '#2b9348', '#55a630', '#80b918', '#aacc00', '#bfd200', '#d4d700', '#dddf00', '#eeef20', '#ffff3f'))
  
  return(res)
}
