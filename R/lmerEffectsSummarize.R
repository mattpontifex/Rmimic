#' lmerEffectsSummarize
#'
#' @description Summary report of univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. Outputted HTML file has additional functionality to expand and collapse results for ease of reading. 
#'
#' @param results list output from the lmerEffects function
#' @param descriptives boolean parameter to control if descriptive statistics are shown before the model result
#' @param modelfit boolean parameter to control if model design and fit statistics are shown
#' @param randomeffects boolean parameter to control if random effects statistics are shown
#' @param fixedeffects boolean parameter to control if fixed effects statistics are shown
#' @param hide list parameter specifying any effects to suppress
#' @param show boolean parameter to control if any results are printed
#' @param outputfile text parameter containing a full file path and file name with html extension. Default is NULL which does not write a file.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#'
#' @importFrom htmlTable concatHtmlTables htmlTable
#' @importFrom purrr list_c
#' @importFrom stats shapiro.test density.default
#' @importFrom stringr str_split
#' @importFrom scales rescale
#' 
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     Rmimic::lmerEffectsSummarize(results, show=TRUE, outputfile=NULL)
#' 
#' @export

lmerEffectsSummarize <- function(results, descriptives=TRUE, modelfit=TRUE, randomeffects=TRUE, fixedeffects=TRUE, style=NULL, hide=NULL, show='html', outputfile=NULL, tag='') {
  
  # since you know within and between, could you use it to create graphs?
  # loop through each effect, bar graphs for between subjects, line graphs for within subjects
  # interactions could 
  
  if (!is.null(style)) {
    if (toupper(style) == toupper("table")) {
      style = "table"
    } else if (toupper(style) == toupper("text")) {
      style = "text"
    }
  } else {
    style = "text"
  }
  boolCSS <- FALSE
  
  if (toupper(show) == toupper("TRUE")) {
    show <- "table"
  } else if (toupper(show) == toupper("html")) {
    show <- "html"
  }
  
  res <- list()
  
  # establish basic style sheet inline
  catout <- ''
  catout <- sprintf('%s%s', catout, Rmimic::RmimicCSSstyle())
  catout <- sprintf('%s<div id="Descriptives"></div>\n', catout)
  res$navtag1 <- catout
  
  standardboardertop <- "border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8;"
  standardboarderbottom <- "border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8;"
  standardpadding <- 'padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; line-height: 1.5em;'
  
  hide <- tryCatch({
    hide <- hide
  }, error = function(e) {
    hide <- c()
  })
  
  studywiseAlpha <- tryCatch({
    studywiseAlpha <- results$studywiseAlpha
  }, error = function(e) {
    studywiseAlpha <- NULL
  })
  if (is.null(studywiseAlpha)) {
    studywiseAlpha <- 0.05
  }
  decimalprecision <- nchar(strsplit(sub('0+$', '', as.character(studywiseAlpha)), ".", fixed = TRUE)[[1]][[2]])
  studywiseAlpharounding <- sort(seq(decimalprecision, (decimalprecision+3), by=1), decreasing=TRUE)
  
  confidenceinterval <- tryCatch({
    confidenceinterval <- results$confidenceinterval
  }, error = function(e) {
    confidenceinterval <- 0.95
  })
  
  if (descriptives) {
    workingdbs <- tryCatch({
      workingdbs <- results$fixedterms
    }, error = function(e) {
      workingdbs <- NULL
    })
    if (!is.null(workingdbs)) {
      # create basic table
      
      if ('descriptives' %in% names(results)) {
        workingdbs <- results$descriptives # if we already have them use them
      } else {
        workingdbs <- Rmimic::lmerEffects_simpledesc(results) # populate them
      }
      
      summarylabel <- ''
      numparticipants <- tryCatch({
        numparticipants <- results$numparticipants
      }, error = function(e) {
        numparticipants <- NULL
      })
      if (!is.null(numparticipants)) {
        if ('meanofparticipants' %in% names(results)) {
          summarylabel <- sprintf('Using %d unique subjects to simulate %d subjects.', numparticipants, results$meanofparticipants)
        } else {
          summarylabel <- sprintf('For %d unique subjects.', numparticipants)
        }
      }
      
      
      catout <- ''
      catout <- sprintf('%s<div style="', catout) 
      catout <- sprintf('%s display: block; height:100%%;', catout)
      catout <- sprintf('%s padding-left:10px;padding-right:10px;padding-top:10px;padding-bottom:10px;', catout)
      catout <- sprintf('%s background-color: #FFFFFF; bgcolor="#FFFFFF"', catout)
      catout <- sprintf('%s" data-quarto-disable-processing="false"; data-quarto-bootstrap="false" width="780" bgcolor="#FFFFFF";>\n', catout)
      
      # add accordian 
      randomsequencetag <- paste(sample(c(letters, LETTERS, 0:100), 32, replace = TRUE), collapse = "")
      catout <- sprintf('%s\n\n<div class="accordion accordion-flush accordion-borderless">\n', catout)
      catout <- sprintf('%s<div class="accordion-item">\n', catout)
      catout <- sprintf('%s<h3 class="accordion-header fw-bold"><button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#idtag%s">\n', catout, randomsequencetag)
    
      # basic table
      catout <- sprintf('%s<table style="display: table;', catout) # start table
      catout <- sprintf('%s table-layout: fixed; width: 100%%; vertical-align: bottom; text-align: center;', catout)
      catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
      catout <- sprintf('%s">\n', catout) # stop table style
      
      # table header
      catout <- sprintf('%s<tr>\n<td style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
    line-height: 1.5;" align="center">', catout) 
      newtextstring <- 'Descriptives <span class="explainmethodtext">(click to expand)</span>'
      catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
      catout <- sprintf('%s<br><span class="secondaryheadings">%s</span>', catout, summarylabel) 
      catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
      catout <- sprintf('%s</table>\n', catout)
      
      # for html output put it in an accordian button
      catout <- sprintf('%s</button></h3>\n', catout) # note that the accordian item is still open at this point
      
      # for html output create accordian body
      catout <- sprintf('%s<div id="idtag%s" class="accordion-collapse collapse">\n', catout, randomsequencetag)
      catout <- sprintf('%s<div class="accordion-body mx-0 px-0">\n<div class="row row-flex px-0 mx-0">\n', catout)
      
      
      colheaders <- c('N', 'Missing', 'Mean', 'Median', 'SD', 'SE', 'Distribution')
      colheadernames <- c('N', 'Missing', 'Mean', 'Median', 'SD', 'SE', 'Distribution')
      numberofcolumns <- length(colheaders)
      
      # basic table
      catout <- sprintf('%s<table style="display: table;', catout) # start table
      catout <- sprintf('%s table-layout: fixed; width: 100%%; vertical-align: bottom; text-align: center;', catout)
      catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
      catout <- sprintf('%s">\n', catout) # stop table style
      
      # headings
      catout <- sprintf('%s<tr style="%s %s %s">\n', catout, standardboardertop, standardboarderbottom, standardpadding) # start row
      for (cC in 1:numberofcolumns) {
        catout <- sprintf('%s<td width="%f%%">', catout, 1/numberofcolumns) # start column
        catout <- sprintf('%s<div class="headinglabeldiv"><span class="basicheadinglabel">%s</span></div>', catout, colheaders[cC]) # start column
        catout <- sprintf('%s</td>\n', catout) # end column
      }
      catout <- sprintf('%s</tr>\n', catout) 
      
      # actual data
      for (cR in 1:nrow(workingdbs)) {
        # set color based upon criteria
        rowcolor <- 'background-color:rgba(255,255,255,0.25);'
        
        # add effect label
        catout <- sprintf('%s<tr style="%s %s %s">\n', catout, rowcolor, standardboardertop, standardpadding) # start row
        catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s"><div class="effectlabeldiv"><span class="effectlabel">%s</span></div></td>\n', catout, (numberofcolumns-1),
                          standardboardertop, workingdbs$Group[cR]) 
        
        # distribution data - creates a stats::density plot in the table 
        catout <- sprintf('%s<td rowspan="2" class="text-center" align="center" style="%s">\n', catout, standardboardertop)
        # color based upon distribution decision
        rectcolor <- '#008934'
        if (workingdbs$DistributionDecision[cR] != "Normal") {
          rectcolor <- '#FFC377'
          #rectcolor <- '#007E8A'
        }
        
        distributioncatout <- tryCatch({
          distributioncatout <- sprintf('<div style="height:100%%; width=100%%;"><svg width="100%%" height="3em">\n')
          datatoplot <- as.numeric(stringr::str_split(workingdbs$DistributionData[cR], ',')[[1]])
          chunkwidth <- 100/length(datatoplot)
          chunkwidthX <- seq(0, 100, length.out=length(datatoplot))
          datatoplot <- (max(datatoplot)-datatoplot)/diff(range(datatoplot))
          datatoplot <- scales::rescale(datatoplot, to=c(0,3.0))
          datatoplotY <- ((datatoplot / 3.0) * 100) * 0.95
          datatoplotHeight <- (100 - datatoplotY) * 3
          for (cDP in 1:length(datatoplot)) {
            distributioncatout <- sprintf('%s<rect width="%f%%" height="%fem" x="%f%%" y="%f%%" fill="%s" style="fill-opacity:0.5;"/>\n',
                                          distributioncatout, chunkwidth, datatoplotHeight[cDP], chunkwidthX[cDP], datatoplotY[cDP], rectcolor)
          }
          distributioncatout <- sprintf('%s</svg></div>\n', distributioncatout)
        }, error = function(e) {
          distributioncatout <- NULL
        })
        if (!is.null(distributioncatout)) {
          catout <- sprintf('%s%s', catout, distributioncatout)
        }
        
        catout <- sprintf('%s</td>\n', catout)
        catout <- sprintf('%s</tr>\n', catout) # end row
        
        # place data
        catout <- sprintf('%s<tr style="%s %s">\n', catout, rowcolor, standardpadding) # start row
        for (cC in 1:(numberofcolumns-1)) {
          catout <- sprintf('%s<td width="%f%%" style="vertical-align: top">', catout, 1/numberofcolumns) # start column
          if (cC < numberofcolumns) {
            catout <- sprintf('%s<div class="demographicdiv"><span class="demographicsstext">%s</span></div>', catout, workingdbs[cR,colheadernames[cC]])
          }
          catout <- sprintf('%s</td>\n', catout) # end column
        }
        catout <- sprintf('%s</tr>\n', catout) # end row
        
      }
      
      catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
      catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
      catout <- sprintf('%s</tr>\n', catout) 
      
      checkindx <- which(workingdbs$DistributionDecision != 'Normal')
      if (length(checkindx) > 0) {
        newtextstring <- "<span style='font-style: italic;'>Note:</span> At least one descriptor was identified as not having normal distribution (see orange distribution)."
        catout <- sprintf('%s<tr><td align="left" colspan="%d"><p class="text-left small">%s</p></td></tr>\n', catout, numberofcolumns, newtextstring)
      }
      
      catout <- sprintf('%s</table></div>\n', catout)
      
      catout <- sprintf('%s</div>\n</div>\n</div>\n</div>\n</div>\n', catout) # closeout accordian
      
      catout <- sprintf('%s<p></p>\n', catout) # end row
      
      
      res$Descriptivesplot <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
    }
  }
  
  catout <- ''
  catout <- sprintf('%s<div id="ModelFit"></div>\n', catout)
  res$navtag2 <- catout
  
  if (modelfit) {
    workingdbs <- tryCatch({
      workingdbs <- results$rsquared
    }, error = function(e) {
      workingdbs <- NULL
    })
    if (!is.null(workingdbs)) {  
      
      # prepare output text
      workingdbs$F2text <- NA
      workingdbs$effects <- round(round(round(round(workingdbs$effects, digits=5), digits=4), digits=3), digits=2)
      workingdbs$ci.lower <- round(round(round(round(workingdbs$ci.lower, digits=5), digits=4), digits=3), digits=2)
      workingdbs$ci.upper <- round(round(round(round(workingdbs$ci.upper, digits=5), digits=4), digits=3), digits=2)
      
      for (cR in 1:nrow(workingdbs)) {
        workingdbs$F2text[cR] <- sprintf('%.2f<br><span style="font-size:80%%; color:#363636;">[%2.0f%% CI: %.2f to %.2f]</span>',
                                         workingdbs$effects[cR], 
                                         floor(results$confidenceinterval*100), workingdbs$ci.lower[cR], workingdbs$ci.upper[cR])
      }
      
      tempworkingdbs <- data.frame(matrix(NA, nrow=1, ncol=3))
      colnames(tempworkingdbs) <- workingdbs$portion
      tempworkingdbs$Fixed[1] <- workingdbs$F2text[1]
      tempworkingdbs$Random[1] <- workingdbs$F2text[2]
      tempworkingdbs$Model[1] <- workingdbs$F2text[3]
      workingdbs <- tempworkingdbs
      
      
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
      
      colheaders <- c('Fixed', "Random", "Model")
      colheadernames <- c('Fixed', "Random", "Model")
      numberofcolumns <- length(colheaders)
      
      # basic table
      catout <- sprintf('%s<table style="display: table;', catout) # start table
      catout <- sprintf('%s table-layout: fixed; width: 100%%; vertical-align: bottom; text-align: center;', catout)
      catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
      catout <- sprintf('%s">\n', catout) # stop table style
      
      # table header
      catout <- sprintf('%s<tr>\n<td colspan="%d" style="padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; 
      line-height: 1.5;" align="center">', catout, numberofcolumns) 
      newtextstring <- 'Summary of Model Fit (R\u200a\u00b2)'
      catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
      catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
      
      
      # headings
      catout <- sprintf('%s<tr style="%s %s %s">\n', catout, standardboardertop, standardboarderbottom, standardpadding) # start row
      for (cC in 1:numberofcolumns) {
        catout <- sprintf('%s<td width="%f%%">', catout, 1/numberofcolumns) # start column
        catout <- sprintf('%s<div class="headinglabeldiv"><span class="basicheadinglabel">%s</span></div>', catout, colheaders[cC]) # start column
        catout <- sprintf('%s</td>\n', catout) # end column
      }
      catout <- sprintf('%s</tr>\n', catout) 
      
      # actual data
      for (cR in 1:nrow(workingdbs)) {
        rowcolor <- 'background-color:rgba(255,255,255,0.25);'
        catout <- sprintf('%s<tr style="height: 5px; %s"></tr>', catout, standardboardertop)
        
        # place data
        catout <- sprintf('%s<tr style="%s %s">\n', catout, rowcolor, standardpadding) # start row
        for (cC in 1:numberofcolumns) {
          catout <- sprintf('%s<td width="%f%%" style="vertical-align: top">', catout, 1/numberofcolumns) # start column
          catout <- sprintf('%s<div class="anovatestresults"><span class="fixedeffectstext">%s</span></div>', catout, workingdbs[cR,colheadernames[cC]])
          catout <- sprintf('%s</td>\n', catout) # end column
        }
        catout <- sprintf('%s</tr>\n', catout) # end row
        
      }
      
      catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
      catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
      catout <- sprintf('%s</tr>\n', catout) 
      
      catout <- sprintf('%s</table></div>\n', catout)
      catout <- sprintf('%s<p></p>\n', catout) # end row
      res$Rsquaredplot <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
      
    }
  }
  
  catout <- ''
  catout <- sprintf('%s<div id="RandomEffects"></div>\n', catout)
  res$navtag3 <- catout
  
  if (randomeffects) {
    workingdbs <- tryCatch({
      workingdbs <- results$randomstats
    }, error = function(e) {
      workingdbs <- NULL
    })
    if (!is.null(workingdbs)) {  
      colnames(workingdbs)[which(colnames(workingdbs) == 'DF')] <- 'df'
      colnames(workingdbs)[which(colnames(workingdbs) == 'LRT')] <- 'Likelihood Ratio'
      colnames(workingdbs)[which(colnames(workingdbs) == 'LogLikelihood')] <- 'log-likelihood'
      colnames(workingdbs)[which(colnames(workingdbs) == 'p.value')] <- 'p'
      
      workingdbs$df <- floor(round(round(round(round(round(workingdbs$df, digits=5), digits=4), digits=3), digits=2), digits=1))
      workingdbs$`log-likelihood` <- round(round(round(round(round(workingdbs$`log-likelihood`, digits=5), digits=4), digits=3), digits=2), digits=1)
      workingdbs$`Likelihood Ratio` <- round(round(round(round(round(workingdbs$`Likelihood Ratio`, digits=5), digits=4), digits=3), digits=2), digits=1)
      workingdbs$p <- sprintf('%.3f', round(round(round(workingdbs$p, digits=5), digits=4), digits=3))
      
      # prepare output text
      for (cR in 1:nrow(workingdbs)) {
        if (workingdbs$p[cR] == "0.000") {
          workingdbs$p[cR] <- '< 0.001'
        }
      }
      
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
      
      if (style == "table") {
        colheaders <- c('', 'df', "Likelihood Ratio", "log-likelihood", "p")
        colheadernames <- c('', 'df', "Likelihood Ratio", "log-likelihood", "p")
      } else {
        colheaders <- c('', '')
        colheadernames <- c('', '')
      }
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
      newtextstring <- 'Summary of Random Effects'
      catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
      catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
      
      # headings
      if (style == "table") {
        catout <- sprintf('%s<tr class="b" style="%s %s %s">\n', catout, standardboardertop, standardboarderbottom, standardpadding) # start row
        catout <- sprintf('%s<td width="%d%%"></td>\n', catout, sigwidth) # flag
        for (cC in 2:numberofcolumns) {
          catout <- sprintf('%s<td width="%f%%">', catout, tempcolwidth) # start column
          catout <- sprintf('%s<div class="headinglabeldiv"><span class="basicheadinglabel">%s</span></div>', catout, colheaders[cC]) # start column
          catout <- sprintf('%s</td>\n', catout) # end column
        }
        catout <- sprintf('%s</tr>\n', catout) 
      }
      
      # actual data
      for (cR in 1:nrow(workingdbs)) {
        # set color based upon criteria
        if (workingdbs$significance[cR]) {
          rowcolor <- 'background-color:rgba(94,201,98,0.25);'
        } else {
          rowcolor <- 'background-color:rgba(255,255,255,0.25);'
        }
        
        # add effect label
        catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
        catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
        catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv"><span class="effectlabel">%s</span></div></td>\n', catout, numberofcolumns-1,
                          standardboardertop, workingdbs$Effect[cR]) 
        catout <- sprintf('%s</tr>\n', catout) # end row
        
        # place data
        if (style == "table") {
          catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
          catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
          for (cC in 2:numberofcolumns) {
            catout <- sprintf('%s<td style="vertical-align: top">', catout) # start column
            catout <- sprintf('%s%s', catout, workingdbs[cR,colheadernames[cC]])
            catout <- sprintf('%s</td>\n', catout) # end column
          }
          catout <- sprintf('%s</tr>\n', catout) # end row
        }
        
        if (!is.na(workingdbs$textoutput[cR])) {
          # add text interpretation
          catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
          catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
          catout <- sprintf('%s<td colspan="%d" align="left">', catout, numberofcolumns-1)
          catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
          catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
          catout <- sprintf('%s<td class="text-center"><div class="anovatestresults"><span class="randomeffectstext">%s</span></div></td>', catout, workingdbs$textoutput[cR])
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
      res$RandomEffectplot <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n\n', catout)
    }
  }
  
  if (fixedeffects) {
    workingdbs <- tryCatch({
      workingdbs <- results$stats
    }, error = function(e) {
      workingdbs <- NULL
    })
    if (!is.null(workingdbs)) {
      
      summarylabel <- sprintf("using %s degrees of freedom approximation", results$df)
      
      workingdbs$DFn <- floor(round(round(round(round(round(workingdbs$DFn, digits=5), digits=4), digits=3), digits=2), digits=1))
      workingdbs$DFd <- floor(round(round(round(round(round(workingdbs$DFd, digits=5), digits=4), digits=3), digits=2), digits=1))
      workingdbs$F.value <- round(round(round(round(round(workingdbs$F.value, digits=5), digits=4), digits=3), digits=2), digits=1)
      workingdbs$p.value <- sprintf('%.3f', round(round(round(workingdbs$p.value, digits=5), digits=4), digits=3))
      workingdbs$fsquared <- sprintf('%.2f', round(round(round(round(workingdbs$fsquared, digits=5), digits=4), digits=3), digits=2))
      workingdbs$fsquared.ci.lower <- sprintf('%.2f', round(round(round(round(workingdbs$fsquared.ci.lower, digits=5), digits=4), digits=3), digits=2))
      workingdbs$fsquared.ci.upper <- sprintf('%.2f', round(round(round(round(workingdbs$fsquared.ci.upper, digits=5), digits=4), digits=3), digits=2))
      
      # prepare output text
      workingdbs$F2text <- NA
      workingdbs$df <- NA
      for (cR in 1:nrow(workingdbs)) {
        workingdbs$F2text[cR] <- sprintf('%s<br><span style="font-size:80%%; color:#363636;">[%2.0f%% CI: %s to %s]</span>',
                                         workingdbs$fsquared[cR], floor(results$confidenceinterval*100), workingdbs$fsquared.ci.lower[cR], workingdbs$fsquared.ci.upper[cR])
        workingdbs$df[cR] <- sprintf('%s, %s', workingdbs$DFn[cR], workingdbs$DFd[cR])
        if (!is.na(workingdbs$textoutput[cR])) {
          workingdbs$textoutput[cR] <- sprintf('<div class="anovatestresults"><span class="fixedeffectstext">%s</span></div>', workingdbs$textoutput[cR])
        }
        if (workingdbs$p.value[cR] == "0.000") {
          workingdbs$p.value[cR] <- '< 0.001'
        }
      }
      
      # check if posthoc
      workingdbs$posthochtml <- NA
      if ('posthoc' %in% names(results)) {
        # get posthoclinkmatrix
        results <- Rmimic::lmerPosthocCorrectionsubprocess(results)
        #results$posthoclinkmatrix
        #results$posthocanovamatrix
        
        # grab posthoc tests
        tempdbs <- Rmimic::lmerEffectsSummarizesubprocess(results)
        workingdbs$posthochtml <- tempdbs$posthochtml
      }
      
      # create navtable
      navtable <- ''
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
      
      
      if (style == "table") {
        colheaders <- c('','df', 'F', 'p', 'cohens f\u200a\u00b2')
        colheadernames <- c('','df', 'F.value', 'p.value', 'F2text')
      } else {
        colheaders <- c('', '')
        colheadernames <- c('', '')
      }
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
      newtextstring <- 'Summary of Fixed Effects'
      catout <- sprintf('%s<p class="text-center"><span class="primaryheadings">%s</span>\n', catout, newtextstring) 
      catout <- sprintf('%s<br><span class="secondaryheadings">%s</span>', catout, summarylabel) 
      catout <- sprintf('%s</p>\n</td>\n</tr>\n', catout) 
      
      # headings
      if (style == "table") {
        catout <- sprintf('%s<tr class="b" style="%s %s %s">\n', catout, standardboardertop, standardboarderbottom, standardpadding) # start row
        catout <- sprintf('%s<td width="%d%%"></td>\n', catout, sigwidth) # flag
        for (cC in 2:numberofcolumns) {
          catout <- sprintf('%s<td width="%f%%">', catout, tempcolwidth) # start column
          catout <- sprintf('%s<div class="headinglabeldiv"><span class="basicheadinglabel">%s</span></div>', catout, colheaders[cC]) # start column
          catout <- sprintf('%s</td>\n', catout) # end column
        }
        catout <- sprintf('%s</tr>\n', catout) 
      }
      
      # actual data
      for (cR in 1:nrow(workingdbs)) {
        # set color based upon criteria
        
        if (workingdbs$significance[cR]) {
          rowcolor <- 'background-color:rgba(94,201,98,0.25);'
        } else {
          rowcolor <- 'background-color:rgba(255,255,255,0.25);'
        }
        
        if (!(workingdbs$Effect[cR] %in% hide)) {
          
          tempnavcol <- 'color: #000;'
          if (!(workingdbs$significance[cR])) {
            tempnavcol <- 'color: #adb5bd;'
          }
          tempnavtable <- sprintf('<li><a class="dropdown-item" style="%s" href="#tnt%s">%s</a></li>\n', tempnavcol, workingdbs$Effect[cR], workingdbs$Effect[cR])
          catout <- sprintf('%s<div id="tnt%s"></div>\n', catout, workingdbs$Effect[cR])
          
          # add effect label
          catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
          catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
          catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv"><span class="effectlabel">%s</span></div></td>\n', catout, numberofcolumns-1,
                            standardboardertop, workingdbs$Effect[cR])
          catout <- sprintf('%s</tr>\n', catout) # end row
          
          # place data
          if (style == "table") {
            catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
            catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
            for (cC in 2:numberofcolumns) {
              catout <- sprintf('%s<td width="%f%%" style="vertical-align: top">', catout, tempcolwidth) # start column
              catout <- sprintf('%s%s', catout, workingdbs[cR,colheadernames[cC]])
              catout <- sprintf('%s</td>\n', catout) # end column
            }
            catout <- sprintf('%s</tr>\n', catout) # end row
          }
          
          if (!is.na(workingdbs$textoutput[cR])) {
            # add text interpretation
            catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
            catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
            catout <- sprintf('%s<td colspan="%d" align="left" style="padding-left: 5px;">', catout, numberofcolumns-1)
            catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
            catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
            catout <- sprintf('%s<td class="text-center">%s</td>', catout, workingdbs$textoutput[cR])
            catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
            catout <- sprintf('%s</td>\n', catout)
            catout <- sprintf('%s</tr>\n', catout) # end row
            
            navtable <- sprintf('%s%s', navtable, tempnavtable)
          }
          if (!is.na(workingdbs$posthochtml[cR])) {
            
            catout <- sprintf('%s<tr style="padding-top: 10px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
            catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
            catout <- sprintf('%s<td colspan="%d" align="left" style="padding-left: 5px;">', catout, numberofcolumns-1)
            catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
            catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
            catout <- sprintf('%s<td class="text-center">%s</td>', catout, workingdbs$posthochtml[cR])
            catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
            catout <- sprintf('%s</td>\n', catout)
            catout <- sprintf('%s</tr>\n', catout) # end row
          }
        } # not hiding
        
      }
      
      catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
      catout <- sprintf('%s<td colspan="%d" class="text-left" align="left"></td>\n', catout, numberofcolumns) 
      catout <- sprintf('%s</tr>\n', catout) 
      
      catout <- sprintf('%s</table>\n', catout)
      catout <- sprintf('%s<p></p>\n', catout) # end row
      
      res$FixedEffectplot <- sprintf('<div class="container">\n<div class="row d-flex justify-content-center">\n<div class="col-md-10 my-2 py-2">%s</div>\n</div>\n</div>\n', catout)
      
    }
  }
  
  
  
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
      textout <- 'lmerEffectsSummarize'
      if (tag != '') {
        textout <- sprintf('%s: %s', textout, tag)
      }
      htmlfilecatout <- sprintf('%s<a class="navbar-brand" href="#">%s</a>\n<button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarNavDropdown" aria-controls="navbarNavDropdown" aria-expanded="false" aria-label="Toggle navigation">\n<span class="navbar-toggler-icon"></span>\n</button>\n\n', htmlfilecatout, textout)
      htmlfilecatout <- sprintf('%s<div class="collapse navbar-collapse" id="navbarNavDropdown">\n<ul class="navbar-nav">\n\n', htmlfilecatout)
      
      if ("Descriptivesplot" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#Descriptives">Descriptives</a></li>\n\n', htmlfilecatout)
      }
      if ("Rsquaredplot" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#ModelFit">Model Fit</a></li>\n\n', htmlfilecatout)
      }
      if ("RandomEffectplot" %in% names(res)) {
        htmlfilecatout <- sprintf('%s<li class="nav-item"><a class="nav-link ps-4" href="#RandomEffects">Random Effects</a></li>\n\n', htmlfilecatout)
      }
      if ("FixedEffectplot" %in% names(res)) {
        if (navtable != '') {
          htmlfilecatout <- sprintf('%s<li class="nav-item dropdown">\n<a class="nav-link ps-4 dropdown-toggle" href="#"  role="button" data-bs-toggle="dropdown" data-bs-auto-close="outside" aria-expanded="false">Fixed Effects</a>\n<ul class="dropdown-menu" aria-labelledby="navbarDropdownMenuLink">\n\n', htmlfilecatout)
          htmlfilecatout <- sprintf('%s%s\n\n', htmlfilecatout, navtable)
          htmlfilecatout <- sprintf('%s</ul>\n\n', htmlfilecatout)
        }
      }
      
      # in future could tag and navigate to posthocs
      #<li class="dropend"><a class="dropdown-item dropdown-toggle" href="#" data-bs-toggle="dropdown" data-bs-auto-close="outside" aria-expanded="false">TitleOfDropDown</a>
      #<ul class="dropdown-menu dropdown-menu-end">
      #  <li><a class="dropdown-item" href="#">Sex</a></li>
      #  <li><a class="dropdown-item" href="#">Group</a></li>
      #</ul>
      #</li>
      
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
    if (show == "table") {
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