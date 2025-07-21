#' lmerRegressionSummarize
#'
#' @description Summary report of linear model with effect size and confidence intervals using a multi-level model from the lme4 function.
#'
#' @param results list output from the lmerModelstats function
#' @param show boolean parameter to control if any results are printed
#' @param outputfile text parameter containing a full file path and file name with html extension. Default is NULL which does not write a file.
#' @param tag text parameter to include additional text in html output file
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom htmlTable concatHtmlTables htmlTable
#' 
#' @examples
#'
#'     altfit <- lmerTest::lmer(Alertness ~ 1 + (1 | PartID), data = Rmimic::alertness)
#'     fit <- lmerTest::lmer(Alertness ~ Group + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerModelstats(fit, altfit)
#'     Rmimic::lmerRegressionSummarize(results, show='html', outputfile=NULL)
#' 
#' @export

lmerRegressionSummarize <- function(results, show='html', outputfile=NULL, tag='') {
 
  boolCSS <- FALSE
  
  res <- list()
  
  # establish basic style sheet inline
  standardboardertop <- "border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8;"
  standardboarderbottom <- "border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8;"
  standardpadding <- 'padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; line-height: 1.5em;'

  
  catout <- ''
  
  # establish basic style sheet inline
  if (!boolCSS) {
    catout <- sprintf('%s%s', catout, Rmimic::RmimicCSSstyle())
    boolCSS <- TRUE
  }
  catout <- sprintf('%s<div style="', catout) 
  catout <- sprintf('%s display: block; height:100%%;', catout)
  catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
  catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
  catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
  
  preheader <- tryCatch({
    preheader <- results$messageout
  }, error = function(e) {
    preheader <- NULL
  })
  if (!is.null(preheader)) {
    catout <- sprintf('%s<p class="text-justify" style="padding-left: 10px; padding-right: 10px;"><span class="preheader">%s</span></p>\n', catout, preheader) # Preheader Information
  }
  catout <- sprintf('%s</div>\n', catout) 
  
  res$preheader <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
  
  res$navtag1 <- sprintf('<div id="RegressionModelFit"></div>\n')
  
  catout <- ''
  
  colheaders <- c('', '')
  colheadernames <- c('', '')
  numberofcolumns <- length(colheaders)
  sigwidth <- 2
  tempcolwidth <- (100-sigwidth)/(length(colheaders) - 1)
  
  catout <- sprintf('%s<div style="', catout) 
  catout <- sprintf('%s display: block; height:100%%;', catout)
  catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
  catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
  catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
  
  # basic table
  catout <- sprintf('%s<table style="display: table;', catout) # start table
  catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
  catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
  catout <- sprintf('%s">\n', catout) # stop table style
  
  # table header
  catout <- sprintf('%s<tr>\n<td colspan="%d" style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
      line-height: 1.5;" align="center">', catout, numberofcolumns) 
  newtextstring <- 'Regression Model Fit'
  catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
  catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
  
  # set color based upon criteria
  if (results$significant) {
    rowcolor <- 'background-color:rgba(94,201,98,0.25);'
  } else {
    rowcolor <- 'background-color:rgba(255,255,255,0.25);'
  }
  
  # add model label
  catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
  catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
  catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv"><span class="effectlabel"><strong>Model:</strong> %s</span></div></td>\n', catout, numberofcolumns-1,
                    standardboardertop, results$fullformula) 
  catout <- sprintf('%s</tr>\n', catout) # end row
  
  # add text 
  catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
  catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
  catout <- sprintf('%s<td colspan="%d" align="left">', catout, numberofcolumns-1)
  catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
  catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
  catout <- sprintf('%s<td class="text-center"><div class="anovatestresults"><span class="fixedeffectstext" style="font-weight: 400; font-size: 0.9em; line-height: 2.00;">%s</span></div></td>', catout, results$text)
  catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
  catout <- sprintf('%s</td>\n', catout)
  catout <- sprintf('%s</tr>\n', catout) # end row
  
  catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
  catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
  catout <- sprintf('%s</tr>\n', catout) 
  
  catout <- sprintf('%s</table></div>\n', catout)
  catout <- sprintf('%s<p></p>\n', catout) # end row
  res$RegressionModelFit <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
  
  res$navtag2 <- sprintf('<div id="RegressionCoefficients"></div>\n')
  
  catout <- ''
  
  # establish basic style sheet inline
  if (!boolCSS) {
    catout <- sprintf('%s%s', catout, Rmimic::RmimicCSSstyle())
    boolCSS <- TRUE
  }
  catout <- sprintf('%s<div style="', catout) 
  catout <- sprintf('%s display: block; height:100%%;', catout)
  catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
  catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
  catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
  
  colheaders <- c('', '')
  colheadernames <- c('', '')
  numberofcolumns <- length(colheaders)
  sigwidth <- 2
  tempcolwidth <- (100-sigwidth)/(length(colheaders) - 1)
  
  # basic table
  catout <- sprintf('%s<table style="display: table;', catout) # start table
  catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
  catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
  catout <- sprintf('%s">\n', catout) # stop table style
  
  # table header
  catout <- sprintf('%s<tr>\n<td colspan="%d" style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
      line-height: 1.5;" align="center">', catout, numberofcolumns) 
  newtextstring <- 'Model Coefficients'
  catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
  catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
  
  # actual data
  for (cR in 1:nrow(results$coefficients)) {
    # set color based upon criteria
    if (results$coefficients$significant[cR]) {
      rowcolor <- 'background-color:rgba(94,201,98,0.25);'
    } else {
      rowcolor <- 'background-color:rgba(255,255,255,0.25);'
    }
    
    # add effect label
    catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
    catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
    catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv"><span class="effectlabel">%s</span></div></td>\n', catout, numberofcolumns-1,
                      standardboardertop, results$coefficients$Variable[cR]) 
    catout <- sprintf('%s</tr>\n', catout) # end row
    
    # place data
    if (!is.na(results$coefficients$text[cR])) {
      # add text interpretation
      catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
      catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
      catout <- sprintf('%s<td colspan="%d" align="left">', catout, numberofcolumns-1)
      catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
      catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
      catout <- sprintf('%s<td class="text-center"><div class="anovatestresults"><span class="randomeffectstext">%s</span></div></td>', catout, results$coefficients$text[cR])
      catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
      catout <- sprintf('%s</td>\n', catout)
      catout <- sprintf('%s</tr>\n', catout) # end row
    }
  }
  
  catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
  catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
  catout <- sprintf('%s</tr>\n', catout) 
  
  catout <- sprintf('%s</table></div>\n', catout)
  catout <- sprintf('%s<p></p>\n', catout) # end row
  res$ModelCoefficients <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
  
  
  
  if (!is.null(results$altfit)) {
  
    res$navtag3 <- sprintf('<div id="RegressionChangeStatistics"></div>\n')
    
    catout <- ''
    
    colheaders <- c('', '')
    colheadernames <- c('', '')
    numberofcolumns <- length(colheaders)
    sigwidth <- 2
    tempcolwidth <- (100-sigwidth)/(length(colheaders) - 1)
    
    catout <- sprintf('%s<div style="', catout) 
    catout <- sprintf('%s display: block; height:100%%;', catout)
    catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
    catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
    catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
    
    # basic table
    catout <- sprintf('%s<table style="display: table;', catout) # start table
    catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
    catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
    catout <- sprintf('%s">\n', catout) # stop table style
    
    # table header
    catout <- sprintf('%s<tr>\n<td colspan="%d" style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
        line-height: 1.5;" align="center">', catout, numberofcolumns) 
    newtextstring <- 'Regression Change Statistics'
    catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
    catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
    
    # set color based upon criteria
    if (results$change$significant) {
      rowcolor <- 'background-color:rgba(94,201,98,0.25);'
    } else {
      rowcolor <- 'background-color:rgba(255,255,255,0.25);'
    }
    
    # add model label
    tempout <- sprintf('<table style="width:100%%">')
    tempout <- sprintf('%s<tr><td style="text-align:right; width:25%%"><span class="effectlabel"><strong>Comparison Model:</strong></span></td><td style="text-align:left;"><span class="effectlabel">%s</span></td></tr>', tempout, results$altformula)
    tempout <- sprintf('%s<tr><td style="text-align:right;"><span class="effectlabel"><strong>Model:</strong></span></td><td style="text-align:left;"><span class="effectlabel">%s</span></td></tr>', tempout, results$fullformula)
    tempout <- sprintf('%s</table>', tempout)
      
    catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
    catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
    catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv">%s</div></td>\n', catout, numberofcolumns-1,
                      standardboardertop, tempout) 
    catout <- sprintf('%s</tr>\n', catout) # end row
    
    # add text 
    catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
    catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
    catout <- sprintf('%s<td colspan="%d" align="left">', catout, numberofcolumns-1)
    catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
    catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
    catout <- sprintf('%s<td class="text-center"><div class="anovatestresults"><span class="fixedeffectstext">%s</span></div></td>', catout, results$changetext)
    catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
    catout <- sprintf('%s</td>\n', catout)
    catout <- sprintf('%s</tr>\n', catout) # end row
    
    catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
    catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
    catout <- sprintf('%s</tr>\n', catout) 
    
    catout <- sprintf('%s</table></div>\n', catout)
    catout <- sprintf('%s<p></p>\n', catout) # end row
    res$RegressionChangeStatistics <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
  }
  
  
  res$navtag4 <- sprintf('<div id="RegressionModelSummary"></div>\n')
  
  catout <- ''
  
  colheaders <- c('', '')
  colheadernames <- c('', '')
  numberofcolumns <- length(colheaders)
  sigwidth <- 2
  tempcolwidth <- (100-sigwidth)/(length(colheaders) - 1)
  
  catout <- sprintf('%s<div style="', catout) 
  catout <- sprintf('%s display: block; height:100%%;', catout)
  catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
  catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
  catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
  
  # basic table
  catout <- sprintf('%s<table style="display: table;', catout) # start table
  catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
  catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
  catout <- sprintf('%s">\n', catout) # stop table style
  
  # table header
  catout <- sprintf('%s<tr>\n<td colspan="%d" style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
      line-height: 1.5;" align="center">', catout, numberofcolumns) 
  newtextstring <- 'Model Summary'
  catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
  catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
  
  # set color based upon criteria
  if (!is.null(results$altfit)) {
    if (results$change$significant) {
      rowcolor <- 'background-color:rgba(94,201,98,0.25);'
    } else {
      rowcolor <- 'background-color:rgba(255,255,255,0.25);'
    }
  } else {
    if (results$significant) {
      rowcolor <- 'background-color:rgba(94,201,98,0.25);'
    } else {
      rowcolor <- 'background-color:rgba(255,255,255,0.25);'
    }    
  }

  # add model label
  catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
  catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
  catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"></td>\n', catout, numberofcolumns-1,
                    standardboardertop) 
  catout <- sprintf('%s</tr>\n', catout) # end row
  
  # add text 
  catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
  catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
  catout <- sprintf('%s<td colspan="%d" align="left">', catout, numberofcolumns-1)
  catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
  catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
  catout <- sprintf('%s<td class="text-left" style=""><div class="anovatestresults"><span class="text-left" style="text-align: left;font-weight: 400; font-size: 0.9em; line-height: 1.75;">%s</span></div></td>', catout, results$modelsummary)
  catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
  catout <- sprintf('%s</td>\n', catout)
  catout <- sprintf('%s</tr>\n', catout) # end row
  
  catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
  catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
  catout <- sprintf('%s</tr>\n', catout) 
  
  catout <- sprintf('%s</table></div>\n', catout)
  catout <- sprintf('%s<p></p>\n', catout) # end row
  res$RegressionModelSummary <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
  
  
  if (length(res) > 0) {
    if (length(res) > 1) {
      catout <- htmlTable::concatHtmlTables(res, headers=rep_len('', length(res)))
    }
    if (!is.null(outputfile)) {
      #writeLines(catout, con = stdout(), sep = "\n", useBytes = FALSE)
      # give it a proper html header
      htmlfilecatout <- sprintf('<!DOCTYPE HTML>\n<html lang="en">\n')
      htmlfilecatout <- sprintf('%s<head>\n<!-- Required meta tags -->\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s<!-- Bootstrap CSS -->\n<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3" crossorigin="anonymous">\n<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s<!-- Masonry JS -->\n<script src="https://cdn.jsdelivr.net/npm/masonry-layout@4.2.2/dist/masonry.pkgd.min.js" integrity="sha384-GNFwBvfVxBkLMJpYMOABq3c+d3KnQxudP/mGPkzpZSTYykLBNsZEnG2D9G/X/+7D" crossorigin="anonymous" async></script>\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s<!-- Jquery JS -->\n<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s<!-- font asesome 6 cdn -->\n<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@fortawesome/fontawesome-free@6.1.0/css/all.css">\n</head>\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s<body>\n\n', htmlfilecatout)
      
      # give page a nav bar
      htmlfilecatout <- sprintf('%s<div class="container">\n<div class="row d-flex justify-content-center my-4">\n<header class="mb-auto">\n<nav class="navbar navbar-expand-lg navbar-dark fixed-top bg-dark" id="navbar-main" >\n<div class="container-fluid">\n\n', htmlfilecatout)
      textout <- 'lmerRegressionSummarize'
      if (tag != '') {
        textout <- sprintf('%s: %s', textout, tag)
      }
      htmlfilecatout <- sprintf('%s<a class="navbar-brand" href="#">%s</a>\n<button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarNavDropdown" aria-controls="navbarNavDropdown" aria-expanded="false" aria-label="Toggle navigation">\n<span class="navbar-toggler-icon"></span>\n</button>\n\n', htmlfilecatout, textout)
      htmlfilecatout <- sprintf('%s<div class="collapse navbar-collapse" id="navbarNavDropdown">\n<ul class="navbar-nav">\n\n', htmlfilecatout)
      
      if ("RegressionModelFit" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#RegressionModelFit">Regression Model Fit</a></li>\n\n', htmlfilecatout)
      }
      if ("ModelCoefficients" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#RegressionCoefficients">Model Coefficients</a></li>\n\n', htmlfilecatout)
      }
      if ("RegressionChangeStatistics" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#RegressionChangeStatistics">Change Statistics</a></li>\n\n', htmlfilecatout)
      }
      if ("RegressionModelSummary" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#RegressionModelSummary">Model Summary</a></li>\n\n', htmlfilecatout)
      }
      
      # end nav bar
      htmlfilecatout <- sprintf('%s</ul>\n</div>\n\n', htmlfilecatout)
      htmlfilecatout <- sprintf('%s</div>\n</nav>\n</header>\n</div>\n</div>\n\n', htmlfilecatout)
      
      # place content
      htmlfilecatout <- sprintf('%s\n%s', htmlfilecatout, catout)
      
      # give it a proper html footer
      htmlfilecatout <- sprintf('%s\n</body>', htmlfilecatout)
      
      fileConn<-file(outputfile)
      writeLines(htmlfilecatout, con = fileConn, sep = "\n", useBytes = FALSE)
      close(fileConn)
    }
    if (show == TRUE) {
      htmlTable::htmlTable(catout)
    }
    if (show == "html") {
      if (!is.null(outputfile)) {
        utils::browseURL(outputfile)
      } else {
        htmlTable::htmlTable(catout)
      }
    }
    #return(writeLines(catout, con = stdout(), sep = "\n", useBytes = FALSE))
  } else {
    return(NULL)
  }
}