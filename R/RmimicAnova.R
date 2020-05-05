#' RmimicAnova
#'
#' @description Compute SPSS style univariate ANOVA with effect size and confidence intervals using the ezANOVA function. Main effects and interactions are automatically decomposed using the specified post-hoc corrections. Subsequent posthoc ANOVAs are outputted as seperate entries in the results for clarity. Note: Data should be in long format.
#'
#' @param data Database containing data
#' @param dependentvariable Dependent Variable label
#' @param subjectid Subject ID label.
#' @param between Between subjects column labels.
#' @param within Within subjects column labels.
#' @param sphericity Parameter to select Sphericity correction method. Default is Greenhouse-Geisser. Other option is Huynh-Feldt.
#' @param feffect Parameter to select which eta square statistic to use for effect size computation. Default is Generalized Eta Squared. Other option is Partial Eta Squared.
#' @param nonparametric Parameter to select use of non parametric t tests. Default is FALSE.
#' @param posthoc Parameter to indicate what post-hoc comparisons should be performed. Default is False Discovery Rate Control. Other options are Bonferroni, Holm-Bonferroni, Scheffe, Sidak, Tukey, or False Discovery Rate Control.
#' @param FDRC Decimal representation of false discocvery rate control. Default is 0.05.
#' @param planned Parameter to specify an effect to show the post-hoc comparisons even if they are not significant.
#' @param suppressposthoc Parameter to specify a posthoc effect to ignore even if it is significant.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#' @param verbosedescriptives Boolean operator for if descriptive statistics should be printed. Default is TRUE.
#' @param posthoclimit Parameter to specify the limit for breaking down interaction terms. Default is 6 indicating a 6 way interaction would not be automatically broken down.
#'
#' @return
#' \item{stats}{ANOVA summary table.}
#' \item{aov}{An aov object corresponding to the requested ANOVA.}
#' \item{posthocttest}{A t test summary table for posthoc decompositions.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 30, 2019
#'
#' @importFrom ez ezANOVA
#' @importFrom doBy summaryBy
#' @importFrom utils packageDate
#' @importFrom pkgcond suppress_conditions
#' 
#'
#' @examples
#'
#'     # Compute ANOVA of PlantGrowth dataset
#'     result <- RmimicAnova(
#'                 data = PlantGrowth,
#'                 dependentvariable = "weight",
#'                 subjectid = NULL,
#'                 between = "group",
#'                 within = NULL)
#'
#'    # Compute ANOVA of ToothGrowth dataset
#'    result <- RmimicAnova(
#'                 data = ToothGrowth,
#'                 dependentvariable = "len",
#'                 subjectid = NULL,
#'                 between = c("dose","supp"),
#'                 within = NULL)
#'
#' @export

RmimicAnova <- function(data, dependentvariable=NULL, subjectid=NULL, between=NULL, within=NULL, sphericity=NULL, feffect=NULL, nonparametric=FALSE, posthoc="False Discovery Rate Control", FDRC=0.05, planned=NULL, suppressposthoc=NULL, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE, verbosedescriptives=TRUE, posthoclimit=6) {

  # revise to incorporate data screening to provide useful error information.
  
  options(contrasts = c("contr.sum", "contr.poly"))
  oldw <- getOption("warn")
  options(warn = -1) # ezANOVA likes to warn about lots of crap

  if (is.null(sphericity)) {
    sphericity <- "Greenhouse-Geisser"
  } else {
    if (toupper(sphericity) == toupper("Huynh-Feldt")) {
      sphericity <- "Huynh-Feldt"
    } else {
      sphericity <- "Greenhouse-Geisser"
    }
  }
  
  if (is.null(feffect)) {
    feffect <- "Generalized Eta Squared"
  } else {
    if (toupper(feffect) == toupper("Partial Eta Squared")) {
      feffect <- "Partial Eta Squared"
    } else if (toupper(feffect) == toupper("Generalized Eta Squared")) {
      feffect <- "Generalized Eta Squared"
    } else {
      feffect <- "Generalized Eta Squared"
    }
  }
  
  if (!is.null(posthoc)) {
    if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
      posthoc <- "False Discovery Rate Control"
    } else if (toupper(posthoc) == toupper("Bonferroni")) {
      posthoc <- "Bonferroni"
    } else if (toupper(posthoc) == toupper("Sidak")) {
      posthoc <- "Sidak"
    } else if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
      posthoc <- "Holm-Bonferroni"
    } else if (posthoc == FALSE) {
      posthoc <- NULL
    } else {
      posthoc <- "False Discovery Rate Control"
    }
  }
  subposthoc <- posthoc
  if (!is.null(posthoc)) {
    if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
      subposthoc <- FALSE
    }
  }
  
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
  
  # Assess what got fed into the function
  betweenvariableL <- length(between)
  withinvariableL <- length(within)
  varlabels <- toupper(colnames(data))
  if (!is.null(subjectid)) {
    data[,subjectid[1]] <- unlist(as.character(data[,subjectid[1]])) # force to characters
    indivparticipant <- sort(unique(as.character(data[,subjectid[1]])))
  } else {
    indivparticipant <- sort(unique(as.character(rownames(data))))
    if ((withinvariableL == 0) & (betweenvariableL > 0)) {
      # no subject ID parameter was provided and the analyis is between subjects
      # assume each is unique
      data$subjectid <- 1:nrow(data)
      subjectid <- "subjectid"
    } else {
      # a within subjects parameter has been specified
      stop("Error in RmimicAnova: A within subjects factor has been requested, but no subjectid parameter was provided to ensure pairwise comparisons. Please provide a column of subject ids.")
    }
  }
  
  completedata <- data[,c(subjectid[1], dependentvariable[1], between, within)]
  completedata <- completedata[complete.cases(completedata),] # remove missing datapoints
  
  # make sure data is collapsed across unnecessary observations
  tempcal <- sprintf("completedata <- doBy::summaryBy(%s ~", dependentvariable[1])
  faclist <- c(subjectid[1], between, within)
  for (cB in 1:length(faclist)) {
    if (cB == 1) {
      tempcal <- sprintf('%s %s', tempcal, faclist[cB])
    } else {
      tempcal <- sprintf('%s + %s', tempcal, faclist[cB])
    }
  }
  tempcal <- sprintf("%s, FUN=c(mean), data=completedata, keep.names=TRUE)", tempcal)
  suppressWarnings(eval(parse(text=tempcal)))
  rm(tempcal, cB)
  
  # Prepare output
  res <- list()
  
  # Compute demographics
  res$demographics <- Rmimic::descriptives(variables=c(dependentvariable[1], between, within), groupvariable = c(between, within), data=completedata, verbose=FALSE)
  
  # Output model
  demoout <- FALSE
  if (verbose == TRUE) {
    temptext <- "Univariate ANOVA Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    # Output model
    outstring <- sprintf("Analysis of %s were conducted using a", dependentvariable[1])
    # Number (Variable: Factors) 
    faclist <- c(between, within)
    for (cB in 1:length(faclist)) {
      factorname <- tolower(faclist[cB])
      factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
      factorsinvolved <- tolower(unique(unlist(as.character(completedata[,faclist[cB]]))))
      outstring <- sprintf("%s %d (%s: %s)", outstring, length(factorsinvolved), factorname, paste(factorsinvolved, collapse = ", "))
      if (cB < length(faclist)) {
        outstring <- sprintf("%s \u00D7", outstring) # seems to work on Mac as well
      }
    }
    outstring <- sprintf('%s univariate', outstring)
    if (!is.null(within)) {
      outstring <- sprintf('%s repeated measures', outstring)
    }
    
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s ANOVA using the ez (Lawrence, %s) and Rmimic (Pontifex, %s) packages in %s.', outstring, strsplit(as.character(utils::packageDate("ez")),"-")[[1]][1],strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1], rvers)
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    # output demographics to console
    if (verbosedescriptives != FALSE) {
      outputdataframe <- res$demographics[,c("CollapsedName", "N", "Mean","Median","SD","SE","Min","Max","Distribution")]
      for (cC in 3:8) {
        outputdataframe[,cC] <- sprintf('%0.1f', round(outputdataframe[,cC], digits=1))
      }
      colnames(outputdataframe)[1] <- 'empty'
      Rmimic::table2console(outputdataframe, sepgap=NULL, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
      rm(outputdataframe)
    }
    
  } # end verbose
  
  # ezANOVA is not presently able to take inputs using dynamic variable names
  funcal <- sprintf('result <- pkgcond::suppress_conditions(ez::ezANOVA(data=completedata,dv=%s,wid=%s,', dependentvariable[1], subjectid[1])
  if (!is.null(between)) {
    tfuncal <- '.('
    for (cB in 1:betweenvariableL) {
      tfuncal <- sprintf('%s%s', tfuncal,between[cB])
      if (cB < betweenvariableL) {
        tfuncal <- sprintf('%s,', tfuncal)
      }
    }
    tfuncal <- sprintf('%s)', tfuncal)
    funcal <- sprintf('%sbetween=%s,', funcal, tfuncal)
    rm(tfuncal)
  } else {
    funcal <- sprintf('%sbetween=NULL,', funcal)
  }
  if (!is.null(within)) {
    tfuncal <- '.('
    for (cB in 1:withinvariableL) {
      tfuncal <- sprintf('%s%s', tfuncal,within[cB])
      if (cB < withinvariableL) {
        tfuncal <- sprintf('%s,', tfuncal)
      }
    }
    rm(cB)
    tfuncal <- sprintf('%s)', tfuncal)
    funcal <- sprintf('%swithin=%s,', funcal, tfuncal)
    rm(tfuncal)
  } else {
    funcal <- sprintf('%s, within=NULL,', funcal)
  }
  funcal <- sprintf('%stype=3,detailed=TRUE,return_aov=TRUE))', funcal)
  
  # Evaluate the text string and tell ezANOVA to shut up, outputs to results
  suppressWarnings(eval(parse(text=funcal)))
  #res$call <- funcal
  rm(funcal)
  
  # obtain effect size estimates and confidence intervals
  result <- Rmimic::ezANOVA2text(result, numparticipants=length(indivparticipant), feffect=feffect, sphericity=sphericity, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha)
  
  res$aov <- result$aov
  res$stats <- result$ANOVA
  
  if (verbose == TRUE) {
    
    # Show ANOVA Table
    if (!is.null(sphericity)) {
      summarylabel <- sprintf("Summary %s corrected ANOVA table", sphericity)
    } else {
      summarylabel <- "Summary ANOVA table"
    }
    cat(sprintf("\n"))
    Rmimic::typewriter(summarylabel, tabs=0, spaces=0, characters=floor(spansize*.9))
    if (toupper(feffect) == toupper("Partial Eta Squared")) {
      summarylabel <- sprintf("With partial eta square based cohens f effect statistics.")
    } else {
      summarylabel <- sprintf("With generalized eta square based cohens f effect statistics.")
    }
    Rmimic::typewriter(summarylabel, tabs=0, spaces=0, characters=floor(spansize*.9))
    
    # prepare output table
    outputdataframe <- data.frame(matrix(NA, nrow=nrow(res$stats), ncol=5))
    # create header labels
    vectnames <- c("Effect", "df", "F", "p")
    if (operatingsystem == "Windows") {
      temptext <- sprintf("f\u00b2 [%2.0f%% CI]", floor(confidenceinterval*100))
    } else {
      temptext <- sprintf("f^2 [%2.0f%% CI]", floor(confidenceinterval*100))
    }
    vectnames <- c(vectnames, temptext)
    colnames(outputdataframe) <- vectnames
    outputdataframe[,1] <- res$stats$Effect
    # manage p value
    for (cR in 1:nrow(outputdataframe)) {
      outputdataframe[cR,3] <- sprintf('%.1f', round(as.numeric(res$stats$F[cR]),digits=1))
      if (outputdataframe[cR,3] == "0.0") {
        outputdataframe[cR,3] <- "< 0.1"
      }
      
      outPvalue <- Rmimic::fuzzyP(res$stats$p[cR])
      if (outPvalue$modifier == "=") {
        pullvalue <- outPvalue$report
      } else {
        pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      if (outPvalue$interpret <= studywiseAlpha) {
        pullvalue <- sprintf('%s%s', pullvalue, "**")
      }
      outputdataframe[cR,4] <- pullvalue
      rm(outPvalue)
    }
    
    # manage df 
    DFn <- round(as.numeric(res$stats$DFn), digits=1)
    DFd <- round(as.numeric(res$stats$DFd), digits=1)
    for (cR in 1:nrow(outputdataframe)) {
      report <- DFn[cR]
      # remove trailing zeros in tens place
      if (substr(report, nchar(report), nchar(report)) == "0") {
        if (substr(report, nchar(report)-1, nchar(report)-1) != ".") {
          report <- substr(report, 1, nchar(report)-1)
        }
      }
      DFn[cR] <- report
      report <- DFd[cR]
      # remove trailing zeros in tens place
      if (substr(report, nchar(report), nchar(report)) == "0") {
        if (substr(report, nchar(report)-1, nchar(report)-1) != ".") {
          report <- substr(report, 1, nchar(report)-1)
        }
      }
      DFd[cR] <- report
      outputdataframe[cR,2] <- sprintf('%s, %s', as.character(DFn[cR]), as.character(DFd[cR])) 
    }
    
    # manage effect size
    for (cR in 1:nrow(outputdataframe)) {
      tempfsq <- sprintf("%.2f", round(as.numeric(res$stats$fsquared[cR]), digits = 2))
      if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
        tempfsq <- "< 0.01"
      }
      tempfsqlw <- sprintf("%.2f", round(as.numeric(res$stats$fsquared.ci.lower[cR]), digits = 2))
      if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
        tempfsqlw <- "0.0"
      }
      tempfsqup <- sprintf("%.2f", round(as.numeric(res$stats$fsquared.ci.upper[cR]), digits = 2))
      if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
        tempfsqup <- "0.0"
      }
      pullvalue <- sprintf('%s [%s, %s]', tempfsq, tempfsqlw, tempfsqup)
      outputdataframe[cR,5] <- pullvalue
      rm(pullvalue)
    }
    
    # Write to console
    Rmimic::table2console(outputdataframe, sepgap=NULL, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)

    # output text writeups
    cat(sprintf("\n\n%s\n", "Summary ANOVA with posthoc decompositions, means, SD, and test statistics"))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  } # end verbose
  
  # populate posthoc decompositions
  posthoclinkmatrix <- data.frame(matrix(NA, nrow=0, ncol=3))
  colnames(posthoclinkmatrix) <- c('p.value', 'Link', 'Row')
  overloadflag <- FALSE
  posthocANOVA <- NULL
  posthocttest <- NULL
  for (currentAnovaLine in 1:nrow(res$stats)) {
    
    effectname <- res$stats$Effect[currentAnovaLine]
    
    # populate a bit more information regarding levels for the ANOVA table
    factorsinvolved <- unlist(strsplit(as.character(res$stats$Effect[currentAnovaLine]),"[:]"))
    factorsinvolvedL <- length(factorsinvolved)
    factorlengthmatrix <- data.frame(matrix(NA,nrow=2, ncol=factorsinvolvedL))
    colnames(factorlengthmatrix) <- factorsinvolved
    for (cB in 1:factorsinvolvedL) {
      templist <- unique(unlist(as.character(completedata[,factorsinvolved[cB]])))
      factorlengthmatrix[1,factorsinvolved[cB]] <- length(templist)
      factorlengthmatrix[2,factorsinvolved[cB]] <- paste(templist, collapse=", ")
    }
    templist <- ''
    for (cB in 1:factorsinvolvedL) {
      templist <- sprintf('%s%s %s (%s)', templist, colnames(factorlengthmatrix)[cB], factorlengthmatrix[1,factorsinvolved[cB]], factorlengthmatrix[2,factorsinvolved[cB]])
      if (cB < factorsinvolvedL) {
        templist <- sprintf('%s : ', templist)
      }
    } 
    res$stats$EffectLevels[currentAnovaLine] <- sprintf('[%s]', templist)
    rm(cB, templist)
    
    # check to see if planned contrast
    forcetrig <- 0
    if (!is.null(planned)) {
      if (planned == TRUE) {
        forcetrig <- 2
      } else {
        for (cXR in 1:length(planned)) {
          if (effectname == planned[cXR]) {
            forcetrig <- 1
          }
        }
        rm(cXR)
      }
    }
    if (!is.null(suppressposthoc)) {
      for (cXR in 1:length(suppressposthoc)) {
        if (effectname == suppressposthoc[cXR]) {
          forcetrig <- -1
        }
      }
      rm(cXR)
    }
    
    # snag p value 
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$stats$p[currentAnovaLine]))
    res$stats$EffectPostHoc[currentAnovaLine] <- 0
    
    
    if (((outPvalue$interpret <= studywiseAlpha) | (forcetrig > 0)) & (!(forcetrig < 0))) {
      # effect was significant or planned
      if (factorsinvolvedL < posthoclimit) {
        
        # subset database for only those factors
        res$stats$EffectPostHoc[currentAnovaLine] <- 1
        workingdatabase <- completedata[,c(subjectid[1], dependentvariable[1], factorsinvolved)]
        
        # make sure data is collapsed across unnecessary observations
        tempcal <- sprintf("workingdatabase <- doBy::summaryBy(%s ~", dependentvariable[1])
        faclist <- c(subjectid[1], factorsinvolved)
        for (cB in 1:length(faclist)) {
          if (cB == 1) {
            tempcal <- sprintf('%s %s', tempcal, faclist[cB])
          } else {
            tempcal <- sprintf('%s + %s', tempcal, faclist[cB])
          }
        }
        tempcal <- sprintf("%s, FUN=c(mean), data=workingdatabase, keep.names=TRUE)", tempcal)
        suppressWarnings(eval(parse(text=tempcal)))
        rm(tempcal, cB)
        
        subbetween <- NULL
        subwithin <- NULL
        for (cB in 1:length(factorsinvolved)) {
          if (factorsinvolved[cB] %in% between) {
            subbetween <- c(subbetween, factorsinvolved[cB])
          }
          if (factorsinvolved[cB] %in% within) {
            subwithin <- c(subwithin, factorsinvolved[cB])
          }
        }
        
        if (factorsinvolvedL == 1) {
          # main effect
          
          # run t-test
          ttestresult <- Rmimic::RmimicTtest(workingdatabase, 
                dependentvariable=dependentvariable[1], 
                subjectid=subjectid[1], 
                between=subbetween, 
                within=subwithin,
                collapse=TRUE, nonparametric=nonparametric, posthoc=subposthoc,
                confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE)
          
          # modify output to indicate subtest
          for (cE in 1:nrow(ttestresult$stats)) {
            ttestresult$stats$EffectNumber[cE] <- currentAnovaLine
            ttestresult$stats$EffectDirection[cE] <- as.character(res$stats$Effect[currentAnovaLine])
          }
          
          # merge output 
          posthocttestL <- 1
          if (is.null(posthocttest)) {
            posthocttest <- ttestresult$stats
          } else {
            posthocttestL <- nrow(posthocttest)+1
            posthocttest <- rbind(posthocttest,ttestresult$stats)
          }
          for (cE in posthocttestL:nrow(posthocttest)) {
            posthoclinkmatrixL <- nrow(posthoclinkmatrix) + 1
            posthoclinkmatrix[posthoclinkmatrixL, 'p.value'] <- posthocttest$p.value[cE]
            posthoclinkmatrix[posthoclinkmatrixL, 'Link'] <- 'posthocttest'
            posthoclinkmatrix[posthoclinkmatrixL, 'Row'] <- cE
          }
          
          rm(subbetween,subwithin, ttestresult)
        } else {
          # interaction
          
          # see what we are working with
          factorlengthmatrix <- data.frame(matrix(NA,nrow=1, ncol=factorsinvolvedL))
          colnames(factorlengthmatrix) <- factorsinvolved
          for (cB in 1:factorsinvolvedL) {
            factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(workingdatabase[,factorsinvolved[cB]]))))
          }
          
          # loop through each factor
          for (currentFactorinvolved in 1:factorsinvolvedL) {
            currentfactorinvolved <- colnames(factorlengthmatrix)[currentFactorinvolved]
            currentfactorlevelsinvolved <- unique(unlist(as.character(workingdatabase[,currentfactorinvolved])))
            otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactorinvolved)]
            
            if (length(otherfactorsinvolved) == 1) {
              decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactorinvolved[1])
            } else {
              decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactorinvolved[1])
            }
            decompdir <- paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 "))
            
            # hold current factor level constant and subset
            for (cD in 1:length(currentfactorlevelsinvolved)) {
              decompconst <- sprintf('[For %s: %s]', currentfactorinvolved[1], currentfactorlevelsinvolved[cD])
              subworkingdatabase <- workingdatabase[which(workingdatabase[,currentfactorinvolved[1]] == currentfactorlevelsinvolved[cD]),]
              
              # make sure data is collapsed across unnecessary observations
              tempcal <- sprintf("subworkingdatabase <- doBy::summaryBy(%s ~", dependentvariable[1])
              faclist <- c(subjectid[1], otherfactorsinvolved)
              for (cE in 1:length(faclist)) {
                if (cE == 1) {
                  tempcal <- sprintf('%s %s', tempcal, faclist[cE])
                } else {
                  tempcal <- sprintf('%s + %s', tempcal, faclist[cE])
                }
              }
              tempcal <- sprintf("%s, FUN=c(mean), data=subworkingdatabase, keep.names=TRUE)", tempcal)
              suppressWarnings(eval(parse(text=tempcal)))
              rm(tempcal, cE)
              
              subbetween <- NULL
              subwithin <- NULL
              for (cE in 1:length(otherfactorsinvolved)) {
                if (otherfactorsinvolved[cE] %in% between) {
                  subbetween <- c(subbetween, otherfactorsinvolved[cE])
                }
                if (otherfactorsinvolved[cE] %in% within) {
                  subwithin <- c(subwithin, otherfactorsinvolved[cE])
                }
              }
              rm(cE)
              
              booltrig <- 0
              # more than one additional variable
              if (length(otherfactorsinvolved) > 1) { 
                booltrig <- 1
              } else {
                # one factor but with more than 2 levels
                otherfactorlevelsinvolved <- unique(unlist(as.character(subworkingdatabase[,otherfactorsinvolved[1]])))
                if (length(otherfactorlevelsinvolved) > 2) { 
                  booltrig <- 1
                }
              }
              
              if (booltrig == 0) {
                # only a single factor with less than 3 levels
                # only a t-test is required
                
                # run t-test
                ttestresult <- Rmimic::RmimicTtest(subworkingdatabase, 
                       dependentvariable=dependentvariable[1], 
                       subjectid=subjectid[1], 
                       between=subbetween, 
                       within=subwithin,
                       collapse=TRUE, nonparametric=nonparametric, posthoc=subposthoc,
                       confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE)

                # modify output to indicate subtest
                for (cE in 1:nrow(ttestresult$stats)) {
                  ttestresult$stats$EffectNumber[cE] <- currentAnovaLine
                  ttestresult$stats$EffectDirection[cE] <- sprintf('%s %s', decompconst, decompdir)
                }
                
                # merge output 
                posthocttestL <- 1
                if (is.null(posthocttest)) {
                  posthocttest <- ttestresult$stats
                } else {
                  posthocttestL <- nrow(posthocttest) + 1
                  posthocttest <- rbind(posthocttest,ttestresult$stats)
                }
                for (cE in posthocttestL:nrow(posthocttest)) {
                  posthoclinkmatrixL <- nrow(posthoclinkmatrix) + 1
                  posthoclinkmatrix[posthoclinkmatrixL, 'p.value'] <- posthocttest$p.value[cE]
                  posthoclinkmatrix[posthoclinkmatrixL, 'Link'] <- 'posthocttest'
                  posthoclinkmatrix[posthoclinkmatrixL, 'Row'] <- cE
                }
                
                rm(subbetween,subwithin, ttestresult)
                
              } else {
                # need to run full anova again
                
                posthocplanned <- NULL
                if (forcetrig > 0) {
                  posthocplanned <- TRUE
                }
                subresult <- Rmimic::RmimicAnova(subworkingdatabase,
                       dependentvariable=dependentvariable[1], 
                       subjectid=subjectid[1], 
                       between=subbetween, 
                       within=subwithin,
                       nonparametric=nonparametric, posthoc=subposthoc,
                       sphericity=sphericity, planned=posthocplanned,
                       verbose=FALSE)
                
                
                # send to output as seperate entry
                res$tempoutputlabel <- subresult
                tempeffectname <- unlist(strsplit(as.character(res$stats$Effect[currentAnovaLine]), split=":"))
                tempeffectname <- paste(tempeffectname, collapse="By")
                if ('posthoclinkmatrix' %in% names(res$tempoutputlabel)) {
                  if (nrow(res$tempoutputlabel$posthoclinkmatrix) > 0) {
                    for (cE in 1:nrow(res$tempoutputlabel$posthoclinkmatrix)) {
                      posthoclinkmatrixL <- nrow(posthoclinkmatrix) + 1
                      posthoclinkmatrix[posthoclinkmatrixL, 'p.value']  <- res$tempoutputlabel$posthoclinkmatrix[cE, 'p.value']
                      posthoclinkmatrix[posthoclinkmatrixL, 'Link']  <- sprintf('PosthocANOVA_%s_%s_%s$%s', tempeffectname, currentfactorlevelsinvolved[cD], paste(otherfactorsinvolved, collapse="By"), res$tempoutputlabel$posthoclinkmatrix[cE, 'Link'])
                      posthoclinkmatrix[posthoclinkmatrixL, 'Row']  <- res$tempoutputlabel$posthoclinkmatrix[cE, 'Row']
                    }
                  }
                }
                names(res)[which(names(res) == 'tempoutputlabel')] <- sprintf('PosthocANOVA_%s_%s_%s', tempeffectname, currentfactorlevelsinvolved[cD], paste(otherfactorsinvolved, collapse="By"))
                
              }
              rm (booltrig)
            }
            rm(cD)
          }
        }
      } else {
        overloadflag <- TRUE
      } # less than 4 factors involved
    } # effect is significant
  } # for each effect
  res$posthocttest <- posthocttest
  res$posthoclinkmatrix <- posthoclinkmatrix
  
  # Perform Posthoc Adjustments
  if (!is.null(posthoc)) {
    #"Bonferroni", "Sidak","Holm-Bonferroni" are addressed within the t-test
    
    if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
      # Glickman, M. E., Rao, S. R., Schultz, M. R. (2014). False discovery rate control is a recommended alternative to Bonferroni-type adjustments in health studies. Journal of Clinical Epidemiology, 67, 850-857.
      
      if ((!is.numeric(FDRC)) | ((FDRC < 0) | (FDRC > 0.99))) {
        FDRC <- 0.05
      }
      if (nrow(res$posthoclinkmatrix) > 0) {
        temp <- posthoclinkmatrix[order(posthoclinkmatrix$p.value),]
        ncomp <- nrow(temp)
        # Loop through P values
        for (rank in 1:nrow(temp)) {
          outPvalue <- Rmimic::fuzzyP(as.numeric(temp$p.value[rank]))
          if (outPvalue$interpret <= studywiseAlpha) {
            temppval <- (FDRC*(rank/ncomp))
            if (as.numeric(temp$p.value[rank]) > as.numeric(temppval)) {
              # P value is no longer considered significant
              criticalphrase <- sprintf("(Benjamini-Hochberg critical alpha = %.3f)", temppval)
              
              templocation <- temp$Link[rank]
              temprow <- temp$Row[rank]
              funcal <- sprintf('tempinterlocation <- res$%s$interpretation[%d]',templocation,temprow)
              suppressWarnings(eval(parse(text=funcal)))
              tempinterlocation <- sprintf("%s However, that difference did not remain significant following false discovery rate control %s.", tempinterlocation, criticalphrase)
              funcal <- sprintf('res$%s$interpretation[%d] <- tempinterlocation',templocation,temprow)
              suppressWarnings(eval(parse(text=funcal)))
              rm(templocation, temprow, funcal, tempinterlocation)
              
            }
          }
        }
        rm(ncomp, temp, rank, outPvalue, temppval, criticalphrase)
      }
    }
  }
  
  if (verbose == TRUE) {
    if (overloadflag == TRUE) {
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      outtext <- sprintf("An interaction exceeding %d variables was detected. To conserve computational resouces, interactions exceeding %d variables should be decomposed in a stepwise fashion manually.", (posthoclimit-1), (posthoclimit-1))
      Rmimic::typewriter(outtext, tabs=0, spaces=0, characters=floor(spansize*.9))
      rm(outtext)
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    }
    
    shortposthoclimit <- (posthoclimit-3)
    if (shortposthoclimit < 3) {
      shortposthoclimit <- 3
    }
    posthoc2text(res, studywiseAlpha=studywiseAlpha, spansize=spansize, currentlevelout=0, posthoclimit=shortposthoclimit, planned=planned, suppressposthoc=suppressposthoc)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  options(warn = oldw) # turn warnings back to original settings
  return(res)
}
  