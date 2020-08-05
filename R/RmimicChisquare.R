#' RmimicChisquare
#'
#' @description Compute SPSS style results for Chi-square analysis with odds ratios and confidence intervals. For samples less than 1000, Fishers exact test statistic is used.
#'
#' @param data Data frame containing the variables of interest.
#' @param x Variable name for the predictor. If x is specified, it overrides the variables input.
#' @param y Variable name for the outcome. If y is specified, it overrides the variables input.
#' @param variables Variable name or list of variables to compute descriptives for. Data should be in the format column 1 is X (predictor), column 2 is Y (outcome).
#' @param posthoc Parameter to indicate what post-hoc comparisons should be performed. Default is False Discovery Rate Control. Other options are Bonferroni, Holm-Bonferroni, Sidak, or False Discovery Rate Control.
#' @param FDRC Decimal representation of false discocvery rate control. Default is 0.05.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param planned Boolean operator for if the posthoc tests should be outputted regardless of the test statistic. Default is FALSE.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#'
#' @return
#' \item{stats}{Summary of analysis.}
#' \item{tabularoutput}{Chi square table.}
#' \item{table}{Simple frequency table.}
#' \item{posthocttest}{A summary table for posthoc tests.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 7, 2020
#'
#' @importFrom epitools oddsratio.fisher oddsratio.wald
#' @importFrom pkgcond suppress_conditions
#' @importFrom utils packageDate
#'
#' @examples
#' 
#'   # was sex related to the chance of survival on the Titanic
#'   chisquareresult <- RmimicChisquare( 
#'           x='Sex', y='Survived', data=Titanic, 
#'           posthoc="False Discovery Rate Control", 
#'           confidenceinterval=0.95, studywiseAlpha=0.05, 
#'           planned=FALSE, verbose=TRUE)
#'          
#'   # is sex related to the likelihood of having a particular hair color 
#'   chisquareresult <- RmimicChisquare( 
#'           x='Sex', y='Hair',  data=HairEyeColor,
#'           posthoc="False Discovery Rate Control", 
#'           confidenceinterval=0.95, studywiseAlpha=0.05, 
#'           planned=FALSE, verbose=TRUE)
#'          
#'   # direct entry dataset example
#'   tempdata <- data.frame("Age"="8","Pet"="Dog","Freq"=282)
#'   tempdata <- rbind(tempdata, data.frame("Age"="30","Pet"="Dog","Freq"=199))
#'   tempdata <- rbind(tempdata, data.frame("Age"="8","Pet"="Cat","Freq"=170))
#'   tempdata <- rbind(tempdata, data.frame("Age"="30","Pet"="Cat","Freq"=240))
#'   chisquareresult <- RmimicChisquare( 
#'           x='Age', y='Pet',  data=tempdata, 
#'           posthoc="False Discovery Rate Control",
#'           confidenceinterval=0.95, studywiseAlpha=0.05, 
#'           planned=FALSE, verbose=TRUE)
#'
#' @export

RmimicChisquare <- function(x=FALSE,y=FALSE,variables=FALSE, data=FALSE, posthoc=TRUE, FDRC=0.05, confidenceinterval=0.95, studywiseAlpha=0.05, planned=FALSE, verbose=TRUE) {  
  
  # debug
  #variables <- c('MotherBMI', 'ChildBMI')
  #data <- workingdata
  #variables=c('carb','cyl')
  #data=mtcars
  #variables=c('treatment','improvement')
  #data=workingdata
  #confidenceinterval=0.95
  #studywiseAlpha=0.05
  #verbose=TRUE
  #planned=TRUE
  
  # Prepare output
  res <- list()
  spansize <- 95
  spancharacter <- "-"
  bigspancharacter <- " - "
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
    bigspancharacter <- " _ "
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
    bigspancharacter <- " - "
  }
  pitchfake <- FALSE

  if (is.table(data)) {
    temp <- as.data.frame(data)
    if (length(which(tolower(names(temp)) == (tolower('freq')))) > 0) {
      # blow it up
      data <- Rmimic::table2frame(data)
    } else {
      Rmimic::typewriter('Alert: Rmimic::RmimicChisquare requires either a data frame or a tabular input.', tabs=0, spaces=2, characters=floor(spansize*.9))
      stop("Rmimic::RmimicChisquare incorrect data input")
    }
  } else {
    if (length(which(tolower(names(data)) == (tolower('freq')))) > 0) {
      # blow it up
      data <- Rmimic::table2frame(data)
    }
  }
  # end table check
  
  if (toupper(posthoc) == toupper("Bonferroni")) {
    posthoc <- "Bonferroni"
  } else if (toupper(posthoc) == toupper("Sidak")) {
    posthoc <- "Sidak"
  } else if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
    posthoc <- "Holm-Bonferroni"
  } else if (posthoc == FALSE) {
    posthoc <- NULL
  } else if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
    posthoc <- "False Discovery Rate Control"
  } else {
    posthoc <- "False Discovery Rate Control"
  }
  
  
  if ((x[1] != FALSE) | (y[1] != FALSE)) {
    if (variables[1] == FALSE) {
      # user specified the x or y call but not the variable call - ideal use
      if ((x[1] != FALSE) & (y[1] != FALSE)) {
        # user specified both x and y
        variables <- c(x[1], y[1])
      } else {
        # uh oh a variable was not entered
        tempvar <- ''
        if (x[1] == FALSE) {
          tempvar <- 'No predictor was specified for the x variable.'
        }
        if (y[1] == FALSE) {
          tempvar <- 'No outcome was specified for the y variable.'
        }
        tempvar <- sprintf('Alert: Rmimic::RmimicChisquare %s.', tempvar)
        Rmimic::typewriter(tempvar, tabs=0, spaces=2, characters=floor(spansize*.9))
        stop("Rmimic::RmimicChisquare incorrect data input")
      }
    } else {
      # user specified the x or y call and the variable call
      if ((x[1] != FALSE) & (y[1] != FALSE)) {
        # both x and y were specified
        variables <- c(x[1], y[1])
      }
    }
  }
  
  if (variables[1] == FALSE) {
    variables <- names(data)
    cDF <- data
  } else {
    cDF <- data.frame(data[,variables])
    names(cDF) <- variables
  }
  cDF <- cDF[which(complete.cases(cDF)),]
  
  # make sure things are factors
  if (!is.factor(cDF[,1])) {
    cDF[,1] <- factor(cDF[,1])
  }
  if (!is.factor(cDF[,1])) {
    cDF[,2] <- factor(cDF[,2])
  }
  
  # verify that the factor levels make sense
  temp1levels <- levels(cDF[,1])
  temp1values <- unique(as.character(unlist(cDF[,1])))
  temp1intersect <- intersect(temp1levels, temp1values)
  newlevels <- c()
  for (cL in 1:length(temp1levels)) {
    if (length(intersect(temp1intersect,temp1levels[cL])) > 0) {
      newlevels <- c(newlevels, temp1levels[cL])
    }
  }
  cDF[,1] <- factor(cDF[,1], levels=newlevels)
  temp1levels <- levels(cDF[,2])
  temp1values <- unique(as.character(unlist(cDF[,2])))
  temp1intersect <- intersect(temp1levels, temp1values)
  newlevels <- c()
  for (cL in 1:length(temp1levels)) {
    if (length(intersect(temp1intersect,temp1levels[cL])) > 0) {
      newlevels <- c(newlevels, temp1levels[cL])
    }
  }
  cDF[,2] <- factor(cDF[,2], levels=newlevels)
  
  
  
  mastertestselect <- sprintf("%s's Exact test", 'Fisher')
  pval <- tryCatch(fisher.test(table(cDF[,1], cDF[,2]), alternative="two.sided")$p.value, error=function(e) NA)
  if (is.na(pval)) {
    mastertestselect <- sprintf("%s's Chi-squared test", 'Pearson')
  }
  # kick to subprocess
  res <- Rmimic::subprocessRmimicChisquare(variables=variables, data=cDF)
  
  if (nrow(cDF) < 1000) {
    # try exact test first, but if there is an error default to Wald
    boolposthoc <- tryCatch({
      results <- pkgcond::suppress_conditions(epitools::oddsratio.fisher(x=cDF[,1], y=cDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
      boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
    )
    if (boolposthoc == FALSE) {
      boolposthoc <- tryCatch({
        results <- pkgcond::suppress_conditions(epitools::oddsratio.wald(x=cDF[,1], y=cDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
        pitchfake <- TRUE
        boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
      )
    }
  } else {
    # Use the Wald test
    boolposthoc <- tryCatch({
      results <- pkgcond::suppress_conditions(epitools::oddsratio.wald(x=cDF[,1], y=cDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
      boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
    )
  }
 
  if (boolposthoc == TRUE) {
    res$table <- results$data # place data table in output
  
    # populate posthoc decompositions
    predictorlevels <- levels(cDF[,1])
    predictorlevelsL <- length(predictorlevels)
    outcomelevels <- levels(cDF[,2])
    outcomelevelsL <- length(outcomelevels)
    
    masterdataframeout <- data.frame(matrix(NA,nrow=0,ncol=10))
    colnames(masterdataframeout) <- c('Reference','Predictor','Outcome','OddsRatio', 'OR.lower.conf.int','OR.upper.conf.int', 'p.value', 'textoutput','overallmodel.report', 'overallmodel.p.value')
    
    # outcomes should binary to make the most sense
    for (cOutcomeLevelsB in 1:(outcomelevelsL)) {
      for (cOutcomeLevelsC in 1:(outcomelevelsL)) {
        if ((cOutcomeLevelsB != cOutcomeLevelsC) & (cOutcomeLevelsB/cOutcomeLevelsC <= 1)) {
      
          # create temporary dataset containing only the binary outcome contrast
          posthoccDF1 <- cDF[which(cDF[,2] == outcomelevels[cOutcomeLevelsB]),]
          posthoccDF2 <- cDF[which(cDF[,2] == outcomelevels[cOutcomeLevelsC]),]
          posthoccDF <- rbind(posthoccDF1,posthoccDF2)
          rm(posthoccDF1,posthoccDF2)
          
          # refactor outcome
          posthoccDF[,2] <- factor(posthoccDF[,2], levels = c(outcomelevels[cOutcomeLevelsB], outcomelevels[cOutcomeLevelsC]))
          
          # rerun on subset
          tempout <- Rmimic::subprocessRmimicChisquare(variables=names(posthoccDF), data=posthoccDF)
          
          # each predictor gets a turn being the reference
          for (cPredictorLevelsB in 1:(predictorlevelsL)) {
            
            dataframeout <- data.frame(matrix(NA,nrow=(predictorlevelsL-1),ncol=10))
            colnames(dataframeout) <- c('Outcome','Reference','Predictor','OddsRatio', 'OR.lower.conf.int','OR.upper.conf.int', 'p.value', 'textoutput', 'overallmodel.report', 'overallmodel.p.value')
            
            subposthoccDF <- posthoccDF
            # refactor predictor
            baselevel <- predictorlevels[cPredictorLevelsB]
            subposthoccDF[,1] <- factor(subposthoccDF[,1], levels = c(baselevel, predictorlevels[which(predictorlevels != baselevel)])) 
            
            if (nrow(subposthoccDF) < 1000) {
              # try exact test first, but if there is an error default to Wald
              boolposthoc <- tryCatch({
                results <- pkgcond::suppress_conditions(epitools::oddsratio.fisher(x=subposthoccDF[,1], y=subposthoccDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
                boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
              )
              if (boolposthoc == FALSE) {
                boolposthoc <- tryCatch({
                  results <- pkgcond::suppress_conditions(epitools::oddsratio.wald(x=subposthoccDF[,1], y=subposthoccDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
                  pitchfake <- TRUE
                  boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
                )
              }
            } else {
              # Use the Wald test
              boolposthoc <- tryCatch({
                results <- pkgcond::suppress_conditions(epitools::oddsratio.wald(x=subposthoccDF[,1], y=subposthoccDF[,2],correction = FALSE,verbose = FALSE,conf.level = confidenceinterval))
                boolposthoc <- TRUE}, error = function(cond){boolposthoc <- FALSE}
              )
            }
            
            if (boolposthoc == TRUE) {
              # cycle through output
              for (cR in 2:(nrow(results$measure))) {
                dataframeout$Outcome[cR-1] <- sprintf('%s vs %s', outcomelevels[cOutcomeLevelsB], outcomelevels[cOutcomeLevelsC])
                dataframeout$Reference[cR-1] <- baselevel
                dataframeout$Predictor[cR-1] <- rownames(results$measure)[cR]
                dataframeout$OddsRatio[cR-1] <- results$measure[cR,1]
                dataframeout$OR.lower.conf.int[cR-1] <- results$measure[cR,2]
                dataframeout$OR.upper.conf.int[cR-1] <- results$measure[cR,3]
                if (nrow(subposthoccDF) < 1000) {
                  dataframeout$p.value[cR-1] <- results$p.value[cR,2]
                  if (is.na(results$p.value[cR,2])) {
                    dataframeout$p.value[cR-1] <- results$p.value[cR,3]
                    if (is.na(results$p.value[cR,3])) {
                      dataframeout$p.value[cR-1] <- results$p.value[cR,1]
                    }
                  }
                } else {
                  dataframeout$p.value[cR-1] <- results$p.value[cR,3]
                  if (is.na(results$p.value[cR,3])) {
                    dataframeout$p.value[cR-1] <- results$p.value[cR,2]
                    if (is.na(results$p.value[cR,2])) {
                      dataframeout$p.value[cR-1] <- results$p.value[cR,1]
                    }
                  }
                }
                outtextstring <- sprintf("Odds Ratio = %.2f", round(dataframeout$OddsRatio[cR-1], digits=2))
                outtextstring <- sprintf("%s [%2.0f%% CI: %.2f to %.2f]", outtextstring, confidenceinterval*100, dataframeout$OR.lower.conf.int[cR-1], dataframeout$OR.upper.conf.int[cR-1])
                outPvalue <- Rmimic::fuzzyP(as.numeric(dataframeout$p.value[cR-1]))
                outtextstring <- sprintf('%s, p %s %s', outtextstring, outPvalue$modifier, outPvalue$report)
                dataframeout$textoutput[cR-1] <- outtextstring
                dataframeout$overallmodel.report[cR-1] <- tempout$stats$textoutput[1]
                dataframeout$overallmodel.p.value[cR-1] <- tempout$stats$p.value[1]
              }
              
              masterdataframeout <- rbind(masterdataframeout, dataframeout)
            } # end posthoc happened
          } # end loop cPredictorLevelsB
        } # prevent duplicates
      } # end loop cOutcomeLevelsC
    } # end loop cOutcomeLevelsB
    res$posthocttest <- masterdataframeout
    
    # Perform Posthoc Adjustments
    if (!is.null(posthoc)) {
      
      # correction should be applied within each outcome chunk, but the nature of the  
      # way comparisons are made allows to just apply it to all contrasts
      if ((toupper(posthoc) == toupper("Bonferroni")) | (toupper(posthoc) == toupper("Sidak"))) {
        critp <- studywiseAlpha
        criticalphrase <- ""
        critdivisor <- length(levels(cDF[,1]))-1
        if (critdivisor < 1) {
          critdivisor <- 1 # shouldnt be possible
        }
        if (toupper(posthoc) == toupper("Bonferroni")) {
          critp <- (studywiseAlpha/critdivisor)
          criticalphrase <- sprintf("(Bonferroni critical alpha = %.4f)", round(critp,4))
        } else if (toupper(posthoc) == toupper("Sidak")) {
          critp <- (1-(1-studywiseAlpha)^(1/critdivisor))
          criticalphrase <- sprintf("(Sidak critical alpha = %.4f)", round(critp,4))
        }
        for (cR in 1:nrow(res$posthocttest)) {
          outPvalue <- Rmimic::fuzzyP(res$posthocttest$p.value[cR])
          if (outPvalue$interpret <= studywiseAlpha) {
            if (outPvalue$interpret > critp) {
              res$posthocttest$textoutput[cR] <- sprintf('%s. However, that effect did not remain significant following correction for multiple comparisons %s', res$posthocttest$textoutput[cR], criticalphrase)
            }
          }
        }
      }
      
      if ((toupper(posthoc) == toupper("Holm-Bonferroni")) | (toupper(posthoc) == toupper("False Discovery Rate Control"))) {
        
        # determine how many actual contrasts are being run and link stats
        temporarycheckframe <- res$posthocttest
        temporarycheckframe$linkrefs <- row.names(temporarycheckframe)
        
        chunks <- unique(res$posthocttest$Outcome)
        masterlinkframe <- data.frame(matrix(NA,nrow=0,ncol=4))
        names(masterlinkframe) <- c('Outcome', 'Contrast', 'p.value', 'Locations')
        
        for (cChunk in 1:length(chunks)) {
          subtemporarycheckframe <- temporarycheckframe[which(temporarycheckframe$Outcome == chunks[cChunk]),] # subset for test chunk
          templevels <- levels(cDF[,1])
          for (cRa in 1:length(templevels)) {
            for (cRb in 1:length(templevels)) {
              if (cRb > cRa) {
                currentlinelist <- c()
                for (cRs in 1:nrow(subtemporarycheckframe)) {
                  if ((subtemporarycheckframe$Reference[cRs] == templevels[cRa]) & (subtemporarycheckframe$Predictor[cRs] == templevels[cRb])) {
                    currentlinelist <- c(currentlinelist, subtemporarycheckframe$linkrefs[cRs])
                  } else if ((subtemporarycheckframe$Reference[cRs] == templevels[cRb]) & (subtemporarycheckframe$Predictor[cRs] == templevels[cRa])) {
                    currentlinelist <- c(currentlinelist, subtemporarycheckframe$linkrefs[cRs])
                  }
                }
                linkframe <- data.frame(matrix(NA,nrow=1,ncol=4))
                names(linkframe) <- c('Outcome', 'Contrast', 'p.value', 'Locations')
                linkframe$Outcome <- chunks[cChunk]
                linkframe$Contrast <- sprintf('%s-%s', templevels[cRa], templevels[cRb])
                if (length(currentlinelist) > 0) {
                  linkframe$p.value <- temporarycheckframe$p.value[as.integer(currentlinelist[1])]
                  linkframe$Locations <- sprintf('%s', paste(currentlinelist, collapse=","))
                }
                masterlinkframe <- rbind(masterlinkframe,linkframe)
              }
            }
          }
        } # end cChunk
        
        
        # Holm-Bonferroni
        if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
          # Holm, S. 1979. A simple sequential rejective multiple test procedure. Scandinavian Journal of Statistics 6:65-70
          
          # should be run within each chunk
          for (cChunk in 1:length(chunks)) {
            submasterlinkframe <- masterlinkframe[which(masterlinkframe$Outcome == chunks[cChunk]),] # subset for test chunk
            submasterlinkframe$Flag <- 0
            
            temp <- data.frame(matrix(NA, nrow=nrow(submasterlinkframe), ncol=2))
            colnames(temp) <- c("p.value", "location")
            for (cR in 1:nrow(submasterlinkframe)) {
              outPvalue <- Rmimic::fuzzyP(submasterlinkframe$p.value[cR])
              if (as.numeric(outPvalue$interpret) <= as.numeric(studywiseAlpha)) {
                temp[cR,1] <- outPvalue$interpret
                temp[cR,2] <- cR
              }
            }
            temp <- temp[order(-temp$p.value),]
            ncomp <- length(which(temp$p.value <= as.numeric(studywiseAlpha)))
            
            # Loop through P values
            rank <- 1
            while (rank <= nrow(temp)) {
              if (!is.na(temp$p.value[rank])) {
                if (as.numeric(temp$p.value[rank]) <= as.numeric(studywiseAlpha)) {
                  temppval <- (studywiseAlpha/(ncomp - rank + 1))
                  if (as.numeric(temp$p.value[rank]) > as.numeric(temppval)) {
                    # P value is no longer considered significant
                    criticalphrase <- sprintf("(Holm-Bonferroni critical alpha = %.4f)", round(temppval, digits=4))
                    submasterlinkframe$Flag[temp$location[rank]] <- sprintf('However, that effect did not remain significant following correction for multiple comparisons %s', criticalphrase)
                  }
                } else {
                  # P value is not significant
                  rank <- rank + 1
                }
              }
              rank <- rank + 1
            }
            
            # add flag to posthoc output
            for (cR in 1:nrow(submasterlinkframe)) {
              if (submasterlinkframe$Flag[cR] != 0) {
                tempvect <- stringr::str_split(submasterlinkframe$Locations[cR], ",")[[1]]
                for (cRc in 1:length(tempvect)) {
                  res$posthocttest$textoutput[as.integer(tempvect[cRc])] <- sprintf('%s. %s', res$posthocttest$textoutput[as.integer(tempvect[cRc])], submasterlinkframe$Flag[cR])
                }
              }
            } # end cR
          } # end cChunk
          
        } # end HB
        
        if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
          # Glickman, M. E., Rao, S. R., Schultz, M. R. (2014). False discovery rate control is a recommended alternative to Bonferroni-type adjustments in health studies. Journal of Clinical Epidemiology, 67, 850-857.
          
          if ((!is.numeric(FDRC)) | ((FDRC < 0) | (FDRC > 0.99))) {
            FDRC <- 0.05
          }
      
          temp <- masterlinkframe[order(masterlinkframe$p.value),]
          ncomp <- nrow(temp)
          # Loop through P values
          for (rank in 1:nrow(temp)) {
            outPvalue <- Rmimic::fuzzyP(as.numeric(temp$p.value[rank]))
            if (outPvalue$interpret <= studywiseAlpha) {
              temppval <- (FDRC*(rank/ncomp))
              if (as.numeric(temp$p.value[rank]) > as.numeric(temppval)) {
                # P value is no longer considered significant
                criticalphrase <- sprintf("(Benjamini-Hochberg critical alpha = %.3f)", temppval)
                tempvect <- stringr::str_split(temp$Locations[rank], ",")[[1]]
                for (cRc in 1:length(tempvect)) {
                  res$posthocttest$textoutput[as.integer(tempvect[cRc])] <- sprintf('%s. However, that effect did not remain significant following false discovery rate control %s', res$posthocttest$textoutput[as.integer(tempvect[cRc])], criticalphrase)
                }
              }
            }
          }
          
        } # FDRC
      } # end HB or FDRC
      
    } # end posthoc corrections
    
  } # posthoc tests got borked
    
  
  if (verbose == TRUE) {
    
    temptext <- "Chi-square Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- "Analysis were conducted using the"
    outstring <- sprintf('%s gmodels (Warnes, Bolker, Lumley, & Johnson, %s)', outstring, strsplit(as.character(utils::packageDate("gmodels")),"-")[[1]][1])
    outstring <- sprintf('%s, epitools (Aragon, %s)', outstring, strsplit(as.character(utils::packageDate("epitools")),"-")[[1]][1])
    outstring <- sprintf('%s, and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    if (pitchfake == FALSE) {
      outstring <- sprintf('%s Test statistics and odds ratios were computed using %s and conditional maximum likelihood estimation', outstring, mastertestselect)
    } else {
      outstring <- sprintf('%s Test statistics and odds ratios were computed using %s and %s unconditional maximum likelihood estimation', outstring, mastertestselect, sprintf("%s's", 'Wald'))
    }
    if ((!is.null(posthoc)) & (posthoc != FALSE)) {
      outstring <- sprintf('%s using the %s approach for post-hoc comparison corrections.', outstring, posthoc)
    }
    
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # show model
    for (cR in 1:length(res$tabularoutput)) {
      cat(sprintf("%s\n",res$tabularoutput[cR]))
    }
    
    outstring <- "Interpretation with Test Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    # report main outcome
    outtext <- ""
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$stats$p.value[1]))
    if (outPvalue$interpret <= studywiseAlpha) {
      outtext <- sprintf("The likelihood of %s was observed to relate to %s", variables[2], variables[1])
    } else {
      outtext <- sprintf("There was no difference in the expected frequency of %s as a function of %s", variables[2], variables[1])
    }
    outtext <- sprintf("%s, %s.", outtext, res$stats$textoutput[1])
    Rmimic::typewriter(outtext, tabs=0, spaces=2, characters=floor(spansize*.9))
    rm(outtext)
    cat(sprintf("\n"))
    
    if (boolposthoc == TRUE) {
    
      if ((outPvalue$interpret <= studywiseAlpha) | (planned == TRUE)) {
        
        tablevel <- 1
        posthocoutcomelevels <- unique(res$posthocttest$Outcome)
        if (length(posthocoutcomelevels) > 1) {
          tablevel <- 2
        }
        for (cPosthocoutcomelevels in 1:length(posthocoutcomelevels)) {
          subposthocDF <- res$posthocttest
          subposthocDF <- subposthocDF[which(subposthocDF$Outcome == posthocoutcomelevels[cPosthocoutcomelevels]),]
          tempvect <- strsplit(posthocoutcomelevels[cPosthocoutcomelevels], " ")[[1]]
          booltrig <- 1
          if (length(posthocoutcomelevels) > 1) {
            booltrig <- 0
            
            Rmimic::typewriter(sprintf('Outcome %s: %s', variables[2], posthocoutcomelevels[cPosthocoutcomelevels]), tabs=1, spaces=0, characters=floor(spansize*.9))
            Rmimic::typewriter(paste(replicate(35, spancharacter), collapse = ""), tabs=1, spaces=0, characters=floor(spansize*.9))
            outtext <- "For clarity, analysis were repeated within each combination of paired outcome levels.\n"
            Rmimic::typewriter(outtext, tabs=1, spaces=0, characters=floor(spansize*.9))
            
            outPvalue2 <- Rmimic::fuzzyP(as.numeric(subposthocDF$overallmodel.p.value[1]))
            if (outPvalue2$interpret <= studywiseAlpha) {
              booltrig <- 1
              outtext <- sprintf("The likelihood of %s (%s, %s) was observed to relate to %s", variables[2], tempvect[1], tempvect[3], variables[1])
            } else {
              outtext <- sprintf("There was no difference in the expected frequency of %s (%s, %s) as function of %s", variables[2], tempvect[1], tempvect[3], variables[1])
            }
            outtext <- sprintf("%s, %s.", outtext, subposthocDF$overallmodel.report[1])
            Rmimic::typewriter(outtext, tabs=1, spaces=2, characters=floor(spansize*.9))
            rm(outtext)
            cat(sprintf("\n"))
          }
          
          if ((booltrig == 1) | (planned == TRUE)) {
              
            posthocpredictorlevels <- unique(subposthocDF$Reference)
            for (cPosthocpredictorlevels in 1:length(posthocpredictorlevels)) {
            
              subsubposthocDF <- subposthocDF[which(subposthocDF$Reference == posthocpredictorlevels[cPosthocpredictorlevels]),]
              
              Rmimic::typewriter(sprintf("%s Breakdown Approach %d %s", paste(replicate(5, spancharacter), collapse = ""), cPosthocpredictorlevels, paste(replicate(5, spancharacter), collapse = "")), tabs=tablevel, spaces=2, characters=floor(spansize*.9))
              decomptext <- sprintf("Post-hoc comparisons were conducted by examining the expected frequency of %s (%s, %s) relative to %s: %s.\n", variables[2], tempvect[1], tempvect[3], variables[1], posthocpredictorlevels[cPosthocpredictorlevels])
              Rmimic::typewriter(decomptext, tabs=tablevel, spaces=2, characters=floor(spansize*.9))
              
              for (cR in 1:nrow(subsubposthocDF)) {
                
                directionalstate <- "an increased"
                if (subsubposthocDF$OddsRatio[cR] <= 1.0) {
                  directionalstate <- "a decreased"
                }
                outPvalue3 <- Rmimic::fuzzyP(as.numeric(subsubposthocDF$p.value[cR]))
                
                tempinsttext <- sprintf('%s %s',posthocpredictorlevels[cPosthocpredictorlevels], variables[1])
                outtext <- sprintf("%s %s", subsubposthocDF$Predictor[cR], variables[1]) 
                if (outPvalue3$interpret <= studywiseAlpha) {
                  outtext <- sprintf("%s was associated with %s likelihood, relative to %s, of %s %s",outtext, directionalstate,tempinsttext,tempvect[3], variables[2])
                } else {
                  outtext <- sprintf("%s made no difference, relative to %s, in the expected frequency of %s %s",outtext,tempinsttext, tempvect[3], variables[2])
                }
                outtext <- sprintf("%s, %s.", outtext, subsubposthocDF$textoutput[cR])
                Rmimic::typewriter(outtext, tabs=tablevel+1, spaces=2, characters=floor(spansize*.9))
                rm(outtext)
                cat(sprintf("\n"))
                
              }
              
              cat(sprintf("\n"))
            } # end loop cPosthocpredictorlevels
          } # booltrig or planned
        } # end loop cPosthocoutcomelevels
      } # if significant or planned
    } # if posthoc test were run
    
    # end output
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  } # end verbose
  
  return(res)
}
  
  