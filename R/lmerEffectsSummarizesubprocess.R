#' lmerEffectsSummarizesubprocess
#'
#' @description A subfunction of lmerEffectsSummarize to place the lmer effects into html.
#'
#' @param res list containing data frames for stats, randomstats, and rsquared values
#' @param mainpasstag text containing a string to tag the results with
#' @param mainfactorlevels numeric indicating the number of test levels
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 29, 2025
#' 
#' @importFrom stats model.frame
#' 
#' @export

lmerEffectsSummarizesubprocess <- function(results, mainpasstag='', mainfactorlevels=0) {
  
  # get posthoclinkmatrix
  results <- Rmimic::lmerPosthocCorrectionsubprocess(results)
  
  tempdbs <- stats::model.frame(results$fit)
  
  standardboardertop <- "border-top-style: solid; border-top-width: 1px; border-top-color: #A8A8A8;"
  standardboarderbottom <- "border-top-style: solid; border-top-width: 1px; border-top-color: #A8A8A8;"
  standardpadding <- 'padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; line-height: 1.5em;'
  
  if ('posthoc' %in% names(results)) {
    
    uniquenames <- names(results$posthoc)
    workingdbs <- results$stats
    workingdbs$posthochtml <- NA
    for (cR in 1:nrow(workingdbs)) {
      
      # see if the posthoc tests were run
      boolchecktestpresent <- FALSE
      if (sprintf('Posthoc_%s', stringr::str_replace_all(workingdbs$Effect[cR], '[:]', '_by_')) %in% uniquenames) {
        boolchecktestpresent <- TRUE
      }
      if (!(boolchecktestpresent)) {
        for (cCheckTestPresent in 1:length(uniquenames)) {
          tempvect <- stringr::str_split(uniquenames[cCheckTestPresent], '_for_')[[1]][1]
          if (sprintf('Posthoc_%s', stringr::str_replace_all(workingdbs$Effect[cR], '[:]', '_by_')) == tempvect) {
            boolchecktestpresent <- TRUE
          }
          if (sprintf('ANOVA_%s', stringr::str_replace_all(workingdbs$Effect[cR], '[:]', '_by_')) == tempvect) {
            boolchecktestpresent <- TRUE
          }
        }
      }
      if (boolchecktestpresent) {
          
        if (workingdbs$factorsinvolved[cR] == 1) {
          # simple t tests
          namecheck <- sprintf('Posthoc_%s', workingdbs$Effect[cR])
          chckindx <- which(results$posthoclinkmatrix$location == namecheck)
          if (length(chckindx) > 0) {
            subposthoclinkmatrix <- results$posthoclinkmatrix[chckindx,]
            
            catout <- '\n\n'
            catout <- sprintf('%s<div style="padding-left: 30px;"><table style="display: table;', catout) # start table
            catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
            catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
            catout <- sprintf('%s">\n', catout) # stop table style
            
            if (nrow(subposthoclinkmatrix) == 1) {
              newtextstring <- '----- Post hoc Comparison -----'
            } else {
              newtextstring <- '----- Post hoc Comparisons -----'
            }
            catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
            catout <- sprintf('%s<td class="text-center" align="center" colspan="2">', catout) # start column
            catout <- sprintf('%s<span class="explainmethodtext">%s</span>', catout, newtextstring) # start column
            catout <- sprintf('%s</td>\n', catout) # end column
            catout <- sprintf('%s</tr>\n', catout) 
            
            for (cContrastOut in 1:nrow(subposthoclinkmatrix)) {
              # set color based upon criteria
              if (subposthoclinkmatrix$significant[cContrastOut]) {
                rowcolor <- 'background-color:rgba(94,201,98,0.25);'
              } else {
                rowcolor <- 'background-color:rgba(255,255,255,0.25);'
              }
              catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
              catout <- sprintf('%s<td class="text-left" align="left" style="%s width: 10px;"></td>\n', catout, rowcolor) 
              catout <- sprintf('%s<td class="text-left" align="left" style=" padding-left: 10px;"><div class="posthoctestresults"><span class="posthocttesteffectstext">%s</span></div></td>\n', catout,
                                subposthoclinkmatrix$textoutput[cContrastOut]) 
              catout <- sprintf('%s</tr>\n', catout) # end row
            }
            
            catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
            catout <- sprintf('%s<td class="text-left" align="left" colspan="2"></td>\n', catout) 
            catout <- sprintf('%s</tr>\n', catout) 
            
            catout <- sprintf('%s</table></div>\n', catout)
            catout <- sprintf('%s\n', catout) # end row
            
            workingdbs$posthochtml[cR] <- catout
          }
        } else {
          # interaction
          factorsinvolved <- unlist(strsplit(as.character(workingdbs$Effect[cR]),"[:]"))
          factorsinvolvedL <- length(factorsinvolved)
          
          # see what we are working with
          factorlengthmatrix <- data.frame(matrix(NA,nrow=1, ncol=factorsinvolvedL))
          colnames(factorlengthmatrix) <- factorsinvolved
          for (cB in 1:factorsinvolvedL) {
            factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(tempdbs[,factorsinvolved[cB]]))))
          }
          
          # loop through each factor
          for (currentFactorinvolved in 1:factorsinvolvedL) {
            currentfactor <- colnames(factorlengthmatrix)[currentFactorinvolved]
            currentfactorlevelsinvolved <- unique(unlist(as.character(tempdbs[,currentfactor])))
            otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactor)]
            
            if (factorsinvolvedL > 1) {
              if (length(otherfactorsinvolved) == 1) {
                decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), otherfactorsinvolved[1], currentfactor[1])
              } else {
                decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), 
                                     paste(otherfactorsinvolved, collapse=sprintf("_by_")), currentfactor[1])
              }
            } else {
              decomptext <- ''
            }
            decompdir <- paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 "))
            
            catout <- '\n\n'
            
            # for html output put it in an accordian button
            # generate a random tag
            randomsequencetag <- paste(sample(c(letters, LETTERS, 0:100), 32, replace = TRUE), collapse = "")
            catout <- sprintf('%s<div class="accordion accordion-flush accordion-borderless">\n', catout)
            catout <- sprintf('%s<div class="accordion-item">\n', catout)
            catout <- sprintf('%s<h3 class="accordion-header fw-bold"><button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#idtag%s">\n', catout, randomsequencetag)

            catout <- sprintf('%s\n<div class="anovadeconstruction">\n<table style="display: table;', catout) # start table
            catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
            catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
            catout <- sprintf('%s">\n', catout) # stop table style
            
            catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
            catout <- sprintf('%s<td class="text-left" align="left" colspan="2">', catout) # start column
            newtextstring <- sprintf('----- Breakdown Approach %d -----', currentFactorinvolved)
            catout <- sprintf('%s<div class="breakdownapproachdiv"><span class="explainmethodtext">%s</span></div>', catout, newtextstring) # start column
            catout <- sprintf('%s</td>\n', catout) # end column
            catout <- sprintf('%s</tr>\n', catout) 
            
            catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
            catout <- sprintf('%s<td class="text-left" align="center" colspan="2">', catout) # start column
            catout <- sprintf('%s<span class="decompdirectiontext">%s</span>', catout, decomptext) # start column
            catout <- sprintf('%s</td>\n', catout) # end column
            catout <- sprintf('%s</tr>\n', catout) 
            
            catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
            catout <- sprintf('%s<td class="text-left" align="left" colspan="2"></td>\n', catout) 
            catout <- sprintf('%s</tr>\n', catout) 
            
            catout <- sprintf('%s</table></div>\n', catout)
            
            # for html output put it in an accordian button
            catout <- sprintf('%s</button></h3>\n', catout) # note that the accordian item is still open at this point
            
            # for html output create accordian body
            catout <- sprintf('%s<div id="idtag%s" class="accordion-collapse collapse">\n', catout, randomsequencetag)
            catout <- sprintf('%s<div class="accordion-body mx-0 px-0">\n<div class="row row-flex px-0 mx-0">\n', catout)
            
            if (!is.na(workingdbs$posthochtml[cR])) {
              workingdbs$posthochtml[cR] <- sprintf('%s\n%s\n', workingdbs$posthochtml[cR], catout) # append text in
            } else {
              workingdbs$posthochtml[cR] <- sprintf('\n%s\n', catout) # place text in
            }
            
            # check if anova is required
            requireanova <- FALSE
            if (factorsinvolvedL > 2) {
              requireanova <- TRUE # 3+ way interaction
            } else {
              # one or two way
              if (length(otherfactorsinvolved) > 0) {
                if (factorlengthmatrix[1,otherfactorsinvolved] > 2) {
                  requireanova <- TRUE # more than 3 levels
                }
              }
            }
            
            if (requireanova) { 
              for (cD in 1:length(currentfactorlevelsinvolved)) {
                # subset data
                #subworkingdatabase <- tempdbs
                #subworkingdatabase <- subworkingdatabase[which(subworkingdatabase[,currentfactor[1]] == currentfactorlevelsinvolved[cD]),]
                
                currentfactorlevelstring <- currentfactorlevelsinvolved[cD]
                # see if factor is just a number
                if (!grepl("\\D", currentfactorlevelstring)) {
                  currentfactorlevelstring <- sprintf('%s%s', currentfactor[1], currentfactorlevelstring)
                }
                decompconst <- sprintf('When %s is %s', currentfactor[1], currentfactorlevelstring)
                
                if (mainfactorlevels == 0) {
                  passtag <- ''
                } else {
                  passtag <- sprintf('&nbsp;&nbsp;%s', mainpasstag)
                }
                
                if (factorsinvolvedL > 1) {
                  fixedformula <- paste(otherfactorsinvolved, collapse=sprintf("*"))
                  if (length(otherfactorsinvolved) == 1) {
                    factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), otherfactorsinvolved[1], currentfactorlevelstring)
                  } else {
                    factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), 
                                         paste(otherfactorsinvolved, collapse=sprintf("_by_")), currentfactorlevelstring)
                  }
                } else {
                  fixedformula <- currentfactor[1]
                  factortag <- currentfactor
                }
                namecheck <- sprintf('ANOVA_%s', factortag)
                if (namecheck %in% uniquenames) {
                  
                  # recursive call
                  textcall <- sprintf("subresults <- results$posthoc$%s", namecheck)
                  eval(parse(text=textcall))
                  subresults <- Rmimic::lmerEffectsSummarizesubprocess(subresults, mainpasstag=sprintf('%s&nbsp;&nbsp;%s', mainpasstag, decompconst), mainfactorlevels=(mainfactorlevels+1))
                 
                  if (!is.null(subresults)) {
                    if (nrow(subresults) > 0) {
                      
                      colheaders <- c('', '')
                      colheadernames <- c('', '')
                      numberofcolumns <- length(colheaders)
                      sigwidth <- 2
                      tempcolwidth <- (100-sigwidth)/(length(colheaders) - 1)
                      
                      # basic table
                      catout <- sprintf('<table style="display: table;') # start table
                      catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
                      catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
                      catout <- sprintf('%s">\n', catout) # stop table style
                      
                      catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
                      catout <- sprintf('%s<td width="%d%%" style=""></td>\n', catout, sigwidth) 
                      catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="padding-left: 5px;"><span class="breakdownconstanttext">%s<span class="explainmethodtext">%s</span></span></td>\n', catout, numberofcolumns-1, decompconst, passtag) 
                      catout <- sprintf('%s</tr>\n', catout) # end row
                      
                      catout <- sprintf('%s</table>\n', catout)
                      catout <- sprintf('%s\n', catout) # end row
                      
                      catout <- sprintf('%s<div style="padding-left: 20px;"><table style="display: table;', catout) # start table
                      catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
                      catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
                      catout <- sprintf('%s">\n', catout) # stop table style
                      
                      
                      passtag <- sprintf('%s&nbsp;&nbsp;%s', mainpasstag, decompconst)
                      
                      # actual data
                      for (cR2 in 1:nrow(subresults)) {
                        # set color based upon criteria
                        
                        if (subresults$significance[cR2]) {
                          rowcolor <- 'background-color:rgba(94,201,98,0.25);'
                        } else {
                          rowcolor <- 'background-color:rgba(255,255,255,0.25);'
                        }
                        
                        # add effect label
                        if (nrow(subresults) > 1) {
                          catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
                          catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
                          catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"><div class="effectlabeldiv"><span class="effectlabel">%s<span class="explainmethodtext">%s</span></span></div></td>\n', catout, numberofcolumns-1,
                                            standardboardertop, subresults$Effect[cR2], passtag)
                          catout <- sprintf('%s</tr>\n', catout) # end row
                        } else {
                          catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
                          catout <- sprintf('%s<td width="%d%%" style="%s %s"></td>\n', catout, sigwidth, rowcolor, standardboardertop) 
                          catout <- sprintf('%s<td colspan="%d" class="text-left" align="left" style="%s padding-left: 5px;"></td>\n', catout, numberofcolumns-1,
                                            standardboardertop)
                          catout <- sprintf('%s</tr>\n', catout) # end row
                        }
                        
                        # place data
                        if (!is.na(subresults$textoutput[cR2])) {
                          # add text interpretation
                          catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
                          catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
                          catout <- sprintf('%s<td colspan="%d" align="left" style="padding-left: 5px;">', catout, numberofcolumns-1)
                          catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
                          catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
                          catout <- sprintf('%s<td class="text-center"><div class="anovatestresults"><span class="fixedeffectstext">%s</span></div></td>', catout, subresults$textoutput[cR2])
                          catout <- sprintf('%s</tr>\n<tr style="height: 10px;"></tr>\n</table>\n', catout)
                          catout <- sprintf('%s</td>\n', catout)
                          catout <- sprintf('%s</tr>\n', catout) # end row
                        }
                        if (!is.na(subresults$posthochtml[cR2])) {
                          catout <- sprintf('%s<tr style="padding-top: 20px; padding-bottom: 20px; padding-left: 5px; padding-right: 5px; border-spacing: 0 50px; line-height: 1.5em;">\n', catout) # start row
                          catout <- sprintf('%s<td style="%s"></td>\n', catout, rowcolor) 
                          catout <- sprintf('%s<td colspan="%d" align="left" style="padding-left: 5px;">', catout, numberofcolumns-1)
                          catout <- sprintf('%s<table align="center"><tr style="width: 100%%; padding-left: 20px;">', catout)
                          catout <- sprintf('%s<td style="width: 10px;"></td>', catout)
                          catout <- sprintf('%s<td class="text-center">%s</td>', catout, subresults$posthochtml[cR2])
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
                      
                      if (!is.na(workingdbs$posthochtml[cR])) {
                        workingdbs$posthochtml[cR] <- sprintf('%s%s', workingdbs$posthochtml[cR], catout) # append text in
                      } else {
                        workingdbs$posthochtml[cR] <- sprintf('%s', catout) # place text in
                      }
                    }
                  }
                } # anova present
              } # cD
            } else {
              # no anova required
              namecheck <- sprintf('Posthoc_%s', factortag)
              chckindx <- which(results$posthoclinkmatrix$location == namecheck)
              if (length(chckindx) > 0) {
                subposthoclinkmatrix <- results$posthoclinkmatrix[chckindx,]
                # make sure nothing weird happened
                if (unique(subposthoclinkmatrix$decomp)[1] == decomptext) {
                  
                  uniqueholds <- unique(subposthoclinkmatrix$hold)
                  for (cUConstant in 1:length(uniqueholds)) {
                    subsubposthoclinkmatrix <- subposthoclinkmatrix[which(subposthoclinkmatrix$hold == uniqueholds[cUConstant]),]
                    decompconst <- sprintf('<span class="breakdownconstanttext">Differences between %s<span class="explainmethodtext">%s</span></span>', stringr::str_replace_all(uniqueholds[cUConstant], '_', ' '), mainpasstag)
                    
                    catout <- '\n\n'
                    catout <- sprintf('%s<div style="padding-left: 30px;"><table style="display: table;', catout) # start table
                    catout <- sprintf('%s width: 100%%; vertical-align: bottom; text-align: center;', catout)
                    catout <- sprintf('%s border-top-style: none; border-bottom-style: none;', catout)
                    catout <- sprintf('%s">\n', catout) # stop table style
                    
                    catout <- sprintf('%s<tr style="%s">\n', catout, standardpadding) # start row
                    catout <- sprintf('%s<td class="text-left" align="left" colspan="2">', catout) # start column
                    catout <- sprintf('%s%s', catout, decompconst) # start column
                    catout <- sprintf('%s</td>\n', catout) # end column
                    catout <- sprintf('%s</tr>\n', catout) 
                    
                    for (cContrastOut in 1:nrow(subsubposthoclinkmatrix)) {
                      # set color based upon criteria
                      if (subsubposthoclinkmatrix$significant[cContrastOut]) {
                        rowcolor <- 'background-color:rgba(94,201,98,0.25);'
                      } else {
                        rowcolor <- 'background-color:rgba(255,255,255,0.25);'
                      }
                      catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboardertop, standardpadding) # start row
                      catout <- sprintf('%s<td class="text-left" align="left" style="%s width: 10px;"></td>\n', catout, rowcolor) 
                      catout <- sprintf('%s<td class="text-left" align="left" style=" padding-left: 10px;"><div class="posthoctestresults"><span class="posthocttesteffectstext">%s</span></div></td>\n', catout,
                                        subsubposthoclinkmatrix$textoutput[cContrastOut]) 
                      catout <- sprintf('%s</tr>\n', catout) # end row
                    }
                    
                    catout <- sprintf('%s<tr style="%s %s">\n', catout, standardboarderbottom, standardpadding) # start row
                    catout <- sprintf('%s<td class="text-left" align="left" colspan="2"></td>\n', catout) 
                    catout <- sprintf('%s</tr>\n', catout) 
                    
                    catout <- sprintf('%s<tr>\n', catout) # start row
                    catout <- sprintf('%s<td colspan="2"></td>\n', catout) 
                    catout <- sprintf('%s</tr>\n', catout) 
                    
                    catout <- sprintf('%s</table></div>\n', catout)
                    catout <- sprintf('%s\n', catout) # end row
                    
                    if (!is.na(workingdbs$posthochtml[cR])) {
                      workingdbs$posthochtml[cR] <- sprintf('%s%s', workingdbs$posthochtml[cR], catout) # append text in
                    } else {
                      workingdbs$posthochtml[cR] <- sprintf('%s', catout) # place text in
                    }
                    
                  } # cUConstant
                } # weird check
              } # posthoc test present
            } # require anova
            
            workingdbs$posthochtml[cR] <- sprintf('%s</div>\n</div>\n</div>\n</div>\n</div>\n', workingdbs$posthochtml[cR]) # closeout accordian
            
            
          } # currentFactorinvolved
          
        } # number of factors
      } # test is present
    } # loop through each name
    return(workingdbs)
  } # posthoc exists
  return(NULL)
}