#' RmimicMLAnova
#'
#' @description Compute SPSS style univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. Main effects and interactions are automatically decomposed using the specified post-hoc corrections. Note: Data should be in long format. As this function will conduct full post-hoc decompositions, the more factors assess the longer processing will take.
#'
#' @param data Database containing data
#' @param dependentvariable Dependent Variable label
#' @param subjectid Subject ID label.
#' @param between Between subjects column labels.
#' @param within Within subjects column labels.
#' @param covariates Non interactive variables to be included in the analysis as covariates.
#' @param randomintercept Parameter to indicate random intercepts.
#' @param randomslope Parameter to indicate random slopes.
#' @param df Parameter to indicate what degrees of freedom approximation should be used. Default is Kenward-Roger. Other option is Shattertwaite.
#' @param posthoc Parameter to indicate what post-hoc comparisons should be performed. Default is False Discovery Rate Control. Other options are Bonferroni, Holm-Bonferroni, Scheffe, Sidak, Tukey, or False Discovery Rate Control.
#' @param FDRC Decimal representation of false discovery rate control. Default is 0.05.
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
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 11, 2020
#'
#' @importFrom utils packageDate
#' @importFrom pkgcond suppress_conditions
#' @importFrom lme4 isSingular
#' @importFrom lmerTest lmer
#' @importFrom emmeans lsmeans
#' @importFrom MuMIn r.squaredGLMM
#'
#'
#' @export

RmimicMLAnova <- function(data, dependentvariable=NULL, subjectid=NULL, between=NULL, within=NULL, covariates=NULL, randomintercept=NULL, randomslope=NULL, df = NULL, posthoc="False Discovery Rate Control", FDRC=0.05, planned=NULL, suppressposthoc=NULL, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE, verbosedescriptives=TRUE, posthoclimit=6) {
  
  #(1|participant) - specifies random intercept, each participant has a line but all lines have the same slope (i.e., the effect is the same for each participant, they just start at different places)
  #(mode|participant) - specifies a random slope for mode for each participant
  
  # NEXT STEPS
  #http://rpsychologist.com/r-guide-longitudinal-lme-lmer
  #https://web.stanford.edu/class/psych253/section/section_8/lmer_examples.html
  #https://www.zoology.ubc.ca/~schluter/R/fit-model/
  
  # debug
  debug <- FALSE
  if (debug == TRUE) {
    data = Rmimic::alertness
    data <- data[which(data$Condition == 'Condition2'),]
    dependentvariable = "Alertness"
    subjectid = "PartID"
    between = "Group"
    #within = c("Condition", "Time")
    within = c("Time")
    covariates = NULL
    randomintercept = c("PartID", "PartID:Condition", "PartID:Time")
    randomintercept = c("PartID")
    randomslope = NULL
    df = NULL
    planned=NULL
    suppressposthoc=NULL
    confidenceinterval=0.95
    studywiseAlpha=0.05
    verbose=TRUE
    verbosedescriptives=TRUE
    posthoclimit=6
    posthoc="False Discovery Rate Control"
    FDRC=0.05
  }
  
  
  #### START FUNCTION
  
  
  options(contrasts = c("contr.sum", "contr.poly"))
  oldw <- getOption("warn")
  options(warn = -1) #
  
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
  
  if (!is.null(df)) {
    if (toupper(df) == toupper("Kenward-Roger")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    }
  } else {
    df = "Kenward-Roger"
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
  
  
  # Assess what got fed into the function
  betweenvariableL <- length(between)
  withinvariableL <- length(within)
  covariateL <- length(covariates)
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
      stop("Error in RmimicMLAnova: A within subjects factor has been requested, but no subjectid parameter was provided to ensure pairwise comparisons. Please provide a column of subject ids.")
    }
  }
  
  # need to handle randomslope and randomintercept to obtain list of factors 
  randomterms <- NULL
  if (!is.null(randomintercept)) {
    for (cB in 1:length(randomintercept)) {
      randomterms <- unique(c(randomterms, unlist(strsplit(as.character(randomintercept[cB]), split=":"))))
    }
  }
  if (!is.null(randomslope)) {
    for (cB in 1:length(randomslope)) {
      randomterms <- unique(c(randomterms, unlist(strsplit(as.character(randomslope[cB]), split="|",fixed=TRUE))))
    }
  }
  fixedterms <- unique(c(between, within))
  allterms <- unique(c(subjectid[1], dependentvariable[1], fixedterms, randomterms, covariates))
  completedata <- data[,allterms] # remove any variables not necessary for the model
  completedata <- completedata[complete.cases(completedata),] # remove missing datapoints
  indivparticipant <- sort(unique(as.character(completedata[,subjectid[1]]))) # recompute
  
  completedata[,subjectid[1]] <- factor(completedata[,subjectid[1]]) # make sure subject is a factor

  # make sure at least two levels of all fixed terms
  erroroutlist <- c()
  for (cB in 1:length(fixedterms)) {
    tempvect <- sort(unique(unlist(as.character(completedata[,fixedterms[cB]]))))
    if (length(tempvect) < 2) {
      erroroutlist <- c(erroroutlist, fixedterms[cB])
    }
  }
  if (length(erroroutlist) > 0) {
    if (length(erroroutlist) > 1) {
      outtext <- sprintf('Alert: Rmimic::RMimicMLAnova the variables %s only contain one level.', paste(erroroutlist, sep=", "))
    } else {
      outtext <- sprintf('Alert: Rmimic::RMimicMLAnova the variable %s only contains a single level.', erroroutlist)
    }
    Rmimic::typewriter(outtext, tabs=0, spaces=2, characters=floor(spansize*.9))
    stop("Rmimic::RMimicMLAnova each factor must have at least two levels")
  }
  
  # Check which random factors can be applied
  booloverfitwarning <- FALSE
  if (!is.null(randomintercept)) {
    randomterms <- NULL
    for (cB in 1:(length(randomintercept))) {
      groupingvariables <- unlist(strsplit(as.character(randomintercept[cB]), split=":"))
      outcal <- 'datamatrix <- expand.grid('
      for (cB2 in 1:length(groupingvariables)) {
        outcal <- sprintf('%s%s = levels(factor(as.character(completedata[,groupingvariables[%d]])))',outcal,groupingvariables[cB2],cB2)
        if (cB2 < length(groupingvariables)) {
          outcal <- sprintf('%s, ', outcal)
        }
      }
      outcal <- sprintf('%s)', outcal)
      eval(parse(text=outcal))
      rm(outcal)
      # number of levels of each grouping factor must be < number of observations
      if (nrow(datamatrix) < nrow(completedata)) {
        randomterms <- c(randomterms, randomintercept[cB])
      }
    }
    randomintercept <- randomterms
    rm(randomterms)
  }
  if (!is.null(randomslope)) {
    randomterms <- NULL
    for (cB in 1:length(randomslope)) {
      temprandomterms <- unlist(strsplit(as.character(tolower(randomslope[cB])), split="|",fixed=TRUE))
      if (temprandomterms[-1] == 1) {
        temprandomterms[1:length(temprandomterms)-1] # remove constant intercepts
      }
      groupingvariables <- c()
      for (cB2 in 1:length(temprandomterms)) {
        groupingvariables <- unique(c(groupingvariables, unlist(strsplit(as.character(temprandomterms[cB2]), split=":"))))
      }
      outcal <- 'datamatrix <- expand.grid('
      for (cB2 in 1:length(groupingvariables)) {
        outcal <- sprintf('%s%s = levels(factor(as.character(completedata[,groupingvariables[%d]])))',outcal,groupingvariables[cB2],cB2)
        if (cB2 < length(groupingvariables)) {
          outcal <- sprintf('%s, ', outcal)
        }
      }
      outcal <- sprintf('%s)', outcal)
      eval(parse(text=outcal))
      rm(outcal)
      # number of levels of each grouping factor must be < number of observations
      if (nrow(datamatrix) < nrow(completedata)) {
        randomterms <- c(randomterms, randomslope[cB])
      }
    }
    randomslope <- randomterms
    rm(randomterms)
  }
  
  # Prepare output
  res <- list()
  
  # Compute demographics
  res$demographics <- Rmimic::descriptives(variables=c(dependentvariable[1]), groupvariable = c(between, within), data=completedata, verbose=FALSE)
  
  # setup random terms
  randomtermsmodel <- c()
  if (!is.null(randomintercept)) {
    for (cB in 1:length(randomintercept)) {
      randomtermsmodel <- c(randomtermsmodel, sprintf("(1 | %s)", randomintercept[cB]))
    }
  }
  if (!is.null(randomslope)) {
    for (cB in 1:length(randomslope)) {
      randomtermsmodel <- c(randomtermsmodel, sprintf("(%s)", randomslope[cB]))
    }
  }
  
  # setup model
  fullmodel <- c()
  factorlength <- 0
  # add non interactive covariates
  if (!is.null(covariates)) {
    fullmodel <- paste(covariates, collapse=" + ")
    factorlength <- factorlength + length(covariates)
  }
  # interactive model of fixed effects
  fullmodel <- paste(c(fullmodel, paste(fixedterms, collapse="*")), collapse=" + ")
  factorlength <- factorlength + length(fixedterms)
  # add random effects
  if (length(randomtermsmodel) > 0 ) {
    fullmodel <- paste(c(fullmodel, paste(randomtermsmodel, collapse=" + ")), collapse=" + ")
    factorlength <- factorlength + length(randomtermsmodel)
  }
  
  # Execute models
  res$model <- sprintf("%s ~ %s",dependentvariable, fullmodel)
  if (length(randomtermsmodel) > 0) {
    
    eval(parse(text=sprintf("fit <- pkgcond::suppress_conditions(lmerTest::lmer(%s, data = completedata))", res$model)))
    # see if the model may be over-fit
    if (lme4::isSingular(fit) == TRUE) {
      booloverfitwarning <- TRUE
    }
    # obtain effect size estimates and confidence intervals
    tempresult <- pkgcond::suppress_conditions(Rmimic::lmer2text(fit, model="ANOVA", df=df, numparticipants=length(indivparticipant), numfactors=factorlength, confidenceinterval=confidenceinterval, data=completedata))
    
    res$fit <- fit
    res$randomstats <- tempresult$RandomEffectsANOVA
    res$stats <- tempresult$ANOVA
    res$Rsquared <- tempresult$Rsquared
  } else {
    Rmimic::typewriter('Alert: Rmimic::RMimicMLAnova insufficient observations within each random effect, reduce the number of random parameters.', tabs=0, spaces=2, characters=floor(spansize*.9))
    stop("Rmimic::RMimicMLAnova number of levels of each random effect exceeds the number of observations")
  }
  
  # Output model
  if (verbose == TRUE) {
    temptext <- "Multi-Level Model Analysis"
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
    outstring <- sprintf('%s multi-level model', outstring) 
    
    if (!is.null(randomintercept)) {
      outstring <- sprintf("%s including the random intercept for", outstring)
      for (cB in 1:length(randomintercept)) {
        temprandomterms <- unlist(strsplit(as.character(tolower(randomintercept[cB])), split=":"))
        for (cB2 in 1:length(temprandomterms)) {
          factorname <- temprandomterms[cB2]
          factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
          outstring <- sprintf("%s %s", outstring, factorname)
          if (length(temprandomterms) > 1) {
            # interaction
            if (cB2 < length(temprandomterms)) {
              outstring <- sprintf("%s \u00D7", outstring) # seems to work on Mac as well
            }
          }
        }
        if (length(randomintercept) > 1) {
          # more than one random term
          if (cB < length(randomintercept)) {
            outstring <- sprintf("%s,", outstring)
          }
          if (cB == (length(randomintercept)-1)) {
            outstring <- sprintf("%s and", outstring)
          }
        }
      }
      if (!is.null(randomslope)) {
        outstring <- sprintf("%s as well as", outstring)
      }
    }
    if (!is.null(randomslope)) {
      if (!is.null(randomintercept)) {
        outstring <- sprintf("%s the random slope for", outstring)
      } else {
        outstring <- sprintf("%s including the random slope for", outstring)
      }
      for (cB in 1:length(randomslope)) {
        temprandomterms <- unlist(strsplit(as.character(tolower(randomslope[cB])), split="|",fixed=TRUE))
        inittemprandomterms <- unlist(strsplit(as.character(tolower(temprandomterms[1])), split=":"))
        if (length(inittemprandomterms) == 1) {
          # main effect in slope
          factorname <- inittemprandomterms[1]
          factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
          outstring <- sprintf("%s %s", outstring, factorname)
        } else {
          # interaction in slope
          for (cB2 in 1:length(inittemprandomterms)) {
            factorname <- inittemprandomterms[cB2]
            factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
            outstring <- sprintf("%s %s", outstring, factorname)
            if (cB2 < length(inittemprandomterms)) {
              outstring <- sprintf("%s \u00D7", outstring) # seems to work on Mac as well
            }
          }
        }
        if (temprandomterms[-1] != 1) {
          outstring <- sprintf("%s for each", outstring)
          inittemprandomterms <- unlist(strsplit(as.character(tolower(temprandomterms[2])), split=":"))
          if (length(inittemprandomterms) == 1) {
            # main effect in slope
            factorname <- inittemprandomterms[1]
            factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
            outstring <- sprintf("%s %s", outstring, factorname)
          } else {
            # interaction in slope
            for (cB2 in 1:length(inittemprandomterms)) {
              factorname <- inittemprandomterms[cB2]
              factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
              outstring <- sprintf("%s %s", outstring, factorname)
              if (cB2 < length(inittemprandomterms)) {
                outstring <- sprintf("%s \u00D7", outstring) # seems to work on Mac as well
              } else {
                outstring <- sprintf("%s interaction", outstring)
              }
            }
          }
        }
        if (length(randomslope) > 1) {
          # more than one random term
          if (cB < length(randomslope)) {
            outstring <- sprintf("%s,", outstring)
          }
          if (cB == (length(randomslope)-1)) {
            outstring <- sprintf("%s and", outstring)
          }
        }
      }
    }
    outstring <- sprintf('%s.', outstring)
    
    if (!is.null(covariates)) {
      outstring <- sprintf('%s Additionally, the model included', outstring)
      for (cB in 1:length(covariates)) {
        factorname <- tolower(covariates[cB])
        factorname <- paste(c(toupper(substr(factorname, 1, 1)), substr(factorname,2,nchar(factorname))), collapse = "")
        factorsinvolved <- tolower(unique(unlist(as.character(completedata[,covariates[cB]]))))
        if (length(factorsinvolved) < 6) {
          outstring <- sprintf("%s %s (%s)", outstring, factorname, paste(factorsinvolved, collapse = ", "))
        } else {
          outstring <- sprintf("%s %s", factorname)
        }
        if (length(covariates) > 1) {
          # more than one random term
          if (cB < length(covariates)) {
            outstring <- sprintf("%s,", outstring)
          }
          if (cB == (length(covariates)-1)) {
            outstring <- sprintf("%s and", outstring)
          }
        } 
      }
      outstring <- sprintf('%s as non-interactive covariates.', outstring)
    }
    
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s Analyses were conducted using the lme4 (Bates et al., %s), lmerTest (Kuznetsova et al., %s), MuMIn (Bartoń, %s), emmeans (Lenth et al., %s), and Rmimic (Pontifex, %s) packages in %s.', outstring, 
                         strsplit(as.character(utils::packageDate("lme4")),"-")[[1]][1],
                         strsplit(as.character(utils::packageDate("lmerTest")),"-")[[1]][1],
                         strsplit(as.character(utils::packageDate("MuMIn")),"-")[[1]][1],
                         strsplit(as.character(utils::packageDate("emmeans")),"-")[[1]][1],
                         strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1], rvers)
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
    outtext <- sprintf('Total unique participants = %d', length(indivparticipant))
    Rmimic::typewriter(outtext, tabs=0, spaces=0, characters=floor(spansize*.9))
    
    # Show ANOVA Table
    summarylabel <- sprintf("Summary of Fixed Effects using %s degrees of freedom approximation", df)
    cat(sprintf("\n"))
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
    
    for (cR in 1:nrow(outputdataframe)) {
      outputdataframe[cR,3] <- sprintf('%.1f', round(as.numeric(res$stats$F.value[cR]),digits=1))
      if (outputdataframe[cR,3] == "0.0") {
        outputdataframe[cR,3] <- "< 0.1"
      }
      
      outPvalue <- Rmimic::fuzzyP(as.numeric(res$stats$p.value[cR]))
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
    
    Rmimic::table2console(outputdataframe, sepgap=NULL, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
    
    tempfsq <- sprintf("%.2f", round(as.numeric(res$Rsquared$FixedEffects[1]), digits = 2))
    if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
      tempfsq <- "< 0.01"
    }
    tempfsqlw <- sprintf("%.2f", round(as.numeric(res$Rsquared$FixedEffects.ci.lower[1]), digits = 2))
    if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
      tempfsqlw <- "0.0"
    }
    tempfsqup <- sprintf("%.2f", round(as.numeric(res$Rsquared$FixedEffects.ci.upper[1]), digits = 2))
    if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
      tempfsqup <- "0.0"
    }
    outtext <- sprintf('Fixed effects r.squared = %s [%2.0f%% CI: %s, %s]', tempfsq, floor(confidenceinterval*100),tempfsqlw,tempfsqup)
    Rmimic::typewriter(outtext, tabs=0, spaces=0, characters=floor(spansize*.9))
    
    
    # Show Random effects ANOVA Table
    summarylabel <- sprintf("Summary of Random Effects")
    cat(sprintf("\n"))
    Rmimic::typewriter(summarylabel, tabs=0, spaces=0, characters=floor(spansize*.9))
    
    outputdataframe <- data.frame(matrix(NA, nrow=nrow(res$randomstats), ncol=5))
    outputdataframe <- res$randomstats[,which(colnames(res$randomstats) != "textoutput")]
    outputdataframe$LogLikelihood <- round(outputdataframe$LogLikelihood, digits=1)
    outputdataframe$LRT <- round(outputdataframe$LRT, digits=1)
    vectnames <- c("Effect", "df", "LogLikelihood", "LRT", "p")
    names(outputdataframe) <-vectnames
    
    for (cR in 1:nrow(outputdataframe)) {
      outPvalue <- Rmimic::fuzzyP(as.numeric(res$randomstats$p.value[cR]))
      if (outPvalue$modifier == "=") {
        pullvalue <- outPvalue$report
      } else {
        pullvalue <- sprintf('%s %s', outPvalue$modifier, outPvalue$report)
      }
      if (outPvalue$interpret <= studywiseAlpha) {
        pullvalue <- sprintf('%s%s', pullvalue, "**")
      }
      outputdataframe[cR,5] <- pullvalue
      rm(outPvalue)
    }
    
    Rmimic::table2console(outputdataframe, sepgap=NULL, spansize=spansize, headers=TRUE, alternate=TRUE, seperators=TRUE)
    
    tempfsq <- sprintf("%.2f", round(as.numeric(res$Rsquared$RandomEffects[1]), digits = 2))
    if ((tempfsq == "= -0.00") | (tempfsq == "= 0.00")) {
      tempfsq <- "< 0.01"
    }
    tempfsqlw <- sprintf("%.2f", round(as.numeric(res$Rsquared$RandomEffects.ci.lower[1]), digits = 2))
    if ((tempfsqlw == "-0.00") | (tempfsqlw == "0.00")) {
      tempfsqlw <- "0.0"
    }
    tempfsqup <- sprintf("%.2f", round(as.numeric(res$Rsquared$RandomEffects.ci.upper[1]), digits = 2))
    if ((tempfsqup == "-0.00") | (tempfsqup == "0.00")) {
      tempfsqup <- "0.0"
    }
    outtext <- sprintf('Random effects r.squared = %s [%2.0f%% CI: %s, %s]', tempfsq, floor(confidenceinterval*100),tempfsqlw,tempfsqup)
    Rmimic::typewriter(outtext, tabs=0, spaces=0, characters=floor(spansize*.9))
    if (booloverfitwarning == TRUE) {
      outtext <- sprintf('Warning: the model is approaching a singular fit, try reducing the number of random effects in the model.')
      Rmimic::typewriter(outtext, tabs=0, spaces=0, characters=floor(spansize*.9))
    }
    
    # output text writeups
    cat(sprintf("\n\n%s\n", "Summary with posthoc decompositions, means, SD, and test statistics"))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  }  
    
  # populate posthoc decompositions
  posthoclinkmatrix <- data.frame(matrix(NA, nrow=0, ncol=3))
  colnames(posthoclinkmatrix) <- c('p.value', 'Link', 'Row')
  overloadflag <- FALSE
  posthocANOVA <- NULL
  posthocttest <- NULL
  res$stats$EffectLevels <- NA
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
    outPvalue <- Rmimic::fuzzyP(as.numeric(res$stats$p.value[currentAnovaLine]))
    res$stats$EffectPostHoc[currentAnovaLine] <- 0
    
    if (((outPvalue$interpret <= studywiseAlpha) | (forcetrig > 0)) & (!(forcetrig < 0))) {
     
      # effect was significant or planned
      if (factorsinvolvedL < posthoclimit) {
        res$stats$EffectPostHoc[currentAnovaLine] <- 1
        
        subbetween <- NULL
        subwithin <- NULL
        for (cB in 1:length(factorsinvolved)) {
          if (factorsinvolved[cB] %in% between) {
            subbetween <- c(subbetween, factorsinvolved[cB])
          }
          if (factorsinvolved[cB] %in% covariates) {
            subbetween <- c(subbetween, factorsinvolved[cB])
          }
          if (factorsinvolved[cB] %in% within) {
            subwithin <- c(subwithin, factorsinvolved[cB])
          }
        }
        
        if (factorsinvolvedL == 1) {
          
          # main effect
          ttestresult <- pkgcond::suppress_conditions(Rmimic::Rmimiclsmeans(fit, completedata, 
                                  dependentvariable=dependentvariable[1], 
                                  subjectid=subjectid[1], 
                                  between=subbetween, 
                                  within=subwithin,
                                  df=df, posthoc=subposthoc,
                                  confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE))
          
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
            factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(completedata[,factorsinvolved[cB]]))))
          }
          
          # loop through each factor
          for (currentFactorinvolved in 1:factorsinvolvedL) {
            currentfactorinvolved <- colnames(factorlengthmatrix)[currentFactorinvolved]
            currentfactorlevelsinvolved <- unique(unlist(as.character(completedata[,currentfactorinvolved])))
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
              subworkingdatabase <- completedata[which(completedata[,currentfactorinvolved[1]] == currentfactorlevelsinvolved[cD]),]
              
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
                #print(sprintf('Nest = %s; Hold = %s', currentfactorinvolved, currentfactorlevelsinvolved[cD]))
                
                ttestresult <- pkgcond::suppress_conditions(Rmimic::Rmimiclsmeans(fit, subworkingdatabase, 
                                             dependentvariable=dependentvariable[1], 
                                             subjectid=subjectid[1], 
                                             between=subbetween, 
                                             within=subwithin,
                                             df=df, posthoc=subposthoc,
                                             nest=currentfactorinvolved,
                                             hold=currentfactorlevelsinvolved[cD],
                                             confidenceinterval=confidenceinterval,studywiseAlpha=studywiseAlpha,verbose=FALSE))
                
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
                
                # Check if there are enough samples for random effects
                if (!is.null(randomintercept)) {
                  randomterms <- NULL
                  for (phcB in 1:(length(randomintercept))) {
                    groupingvariables <- unlist(strsplit(as.character(randomintercept[phcB]), split=":"))
                    outcal <- 'datamatrix <- expand.grid('
                    for (phcB2 in 1:length(groupingvariables)) {
                      outcal <- sprintf('%s%s = levels(factor(as.character(subworkingdatabase[,groupingvariables[%d]])))',outcal,groupingvariables[phcB2],phcB2)
                      if (phcB2 < length(groupingvariables)) {
                        outcal <- sprintf('%s, ', outcal)
                      }
                    }
                    outcal <- sprintf('%s)', outcal)
                    eval(parse(text=outcal))
                    rm(outcal)
                    # number of levels of each grouping factor must be < number of observations
                    if (nrow(datamatrix) < nrow(subworkingdatabase)) {
                      randomterms <- c(randomterms, randomintercept[phcB])
                    }
                  }
                  randomintercept <- randomterms
                  rm(randomterms)
                }
                if (!is.null(randomslope)) {
                  randomterms <- NULL
                  for (phcB in 1:length(randomslope)) {
                    temprandomterms <- unlist(strsplit(as.character(tolower(randomslope[phcB])), split="|",fixed=TRUE))
                    if (temprandomterms[-1] == 1) {
                      temprandomterms[1:length(temprandomterms)-1] # remove constant intercepts
                    }
                    groupingvariables <- c()
                    for (phcB2 in 1:length(temprandomterms)) {
                      groupingvariables <- unique(c(groupingvariables, unlist(strsplit(as.character(temprandomterms[phcB2]), split=":"))))
                    }
                    outcal <- 'datamatrix <- expand.grid('
                    for (phcB2 in 1:length(groupingvariables)) {
                      outcal <- sprintf('%s%s = levels(factor(as.character(subworkingdatabase[,groupingvariables[%d]])))',outcal,groupingvariables[phcB2],phcB2)
                      if (phcB2 < length(groupingvariables)) {
                        outcal <- sprintf('%s, ', outcal)
                      }
                    }
                    outcal <- sprintf('%s)', outcal)
                    eval(parse(text=outcal))
                    rm(outcal)
                    # number of levels of each grouping factor must be < number of observations
                    if (nrow(datamatrix) < nrow(subworkingdatabase)) {
                      randomterms <- c(randomterms, randomslope[phcB])
                    }
                  }
                  randomslope <- randomterms
                  rm(randomterms)
                }
                
                if (is.null(randomintercept) & is.null(randomslope)) {
                  
                  #print(sprintf('Line844 - repeat full model but use RmimicAnova'))
                  # have to run an ANOVA then
                  subresult <- pkgcond::suppress_conditions(Rmimic::RmimicAnova(subworkingdatabase,
                                   dependentvariable=dependentvariable[1], 
                                   subjectid=subjectid[1], 
                                   between=subbetween, 
                                   within=subwithin,
                                   posthoc=subposthoc,
                                   planned=posthocplanned,
                                   confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha,
                                   verbose=FALSE))
                  
                } else {
                  
                  #print(sprintf('Line844 - repeat full model but use RmimicMLAnova'))
                  subresult <- pkgcond::suppress_conditions(Rmimic::RmimicMLAnova(subworkingdatabase,
                                   dependentvariable=dependentvariable[1], 
                                   subjectid=subjectid[1], 
                                   between=subbetween, 
                                   within=subwithin,
                                   covariates=covariates, 
                                   randomintercept=randomintercept, randomslope=randomslope, 
                                   df = df, posthoc=subposthoc, 
                                   planned=posthocplanned, 
                                   confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, 
                                   verbose=FALSE))
                  
                }
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
                
              } # end booltrig
            } # end cD
          } # end currentFactorinvolved
        } # end factorsinvolvedL
      } else {
        overloadflag <- TRUE
      } # end have not exceeded posthoc limit
    } # end significant or planned
  } # end currentAnovaLine
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
        rm(ncomp, temp, rank, outPvalue, temppval)
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
    Rmimic::posthoc2text(res, studywiseAlpha=studywiseAlpha, spansize=spansize, currentlevelout=0, posthoclimit=shortposthoclimit, planned=planned, suppressposthoc=suppressposthoc)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  options(warn = oldw) # turn warnings back to original settings
  return(res)
}
  
  
