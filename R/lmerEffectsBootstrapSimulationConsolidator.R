#' lmerEffectsBootstrapSimulationConsolidator
#'
#' @description Wrapper to perform bootstrap analysis on an lmer object processed through the Rmimic::lmerEffects function. This portion handles the consolidation of the models and summarizes them into a singular model.
#'
#' @param results List containing output from lmerEffects
#' @param average text parameter to indicate how the results should be collapsed. Default is using the median. Other options is mean.
#' @param reporteddata text parameter to indicate if the posthoc text reports should use actual data (actual) or should report the simulated or resampled data means.
#'
#' @return
#' \item{results}{List containing boostrapped output.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 23, 2025
#'
#' @importFrom Rmisc CI
#' @importFrom miscTools colMedians
#' @importFrom stats model.frame
#' @importFrom progressr progressor with_progress
#' 
#' @examples
#' \dontrun{
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
#'                covariates=NULL, planned=c('Group'), posthoccorrection="False Discovery Rate Control", progressbar=TRUE)
#'     results <- Rmimic::lmerEffectsBootstrapSimulationCreation(results, repetitions=999, tmpdir='/tempdirectory/')
#'     results <- Rmimic::lmerEffectsBootstrapSimulationConsolidator(results, average='median', reporteddata='actual')
#'     }
#'
#' @export

lmerEffectsBootstrapSimulationConsolidator <- function(results, average='median', reporteddata='simulated') {
  
  invisible(suppressWarnings(suppressPackageStartupMessages(suppressMessages(library(progressr)))))
  
  # with large numbers of files, we end up with a memory issue
  # one solution is a disk-based approach
  
  # create folder to hold models
  tmpdir <- file.path(results$futuretag$tmpdir, sprintf('modelsummary_%s', paste(sample(c(0:9, letters, LETTERS), 10, replace=TRUE), collapse="")))
  dir.create(tmpdir, showWarnings = FALSE)
  
  # create list of outputs
  find_tables <- function(x, prefix = "") {
    find_tables_results <- character()
    if (is.list(x)) {
      for (name in names(x)) {
        current_name <- if (prefix == "") name else paste0(prefix, "$", name)
        element <- x[[name]]
        if (is.data.frame(element)) {
          find_tables_results <- c(find_tables_results, current_name)
        } else if (is.list(element)) {
          find_tables_results <- c(find_tables_results, find_tables(element, current_name))
        }
      }
    } else if (is.data.frame(x)) {
      find_tables_results <- c(find_tables_results, prefix)
    }
    return(find_tables_results)
  }
  elementsofinterest <- c('stats', 'randomstats', 'rsquared', 'descriptives', 'numparticipants')
  if ('posthoc' %in% names(results)) {
    tempvect <- find_tables(results$posthoc)
    tempvect <- paste0("posthoc$", tempvect)
    elementsofinterest <- c(elementsofinterest, tempvect)
  }
  tablelookupdirectory <- data.frame('Real' = elementsofinterest)
  
  generate_unique_strings <- function(n, string_length = 20, pool = c(0:9, letters, LETTERS)) {
    unique_strings <- character()
    while (length(unique_strings) < n) {
      needed <- n - length(unique_strings)
      new_strings <- replicate(needed, paste(sample(pool, string_length, replace = TRUE), collapse = ""))
      unique_strings <- unique(c(unique_strings, new_strings))
    }
    unique_strings <- paste0(unique_strings, ".RData")
    return(unique_strings)
  }
  tablelookupdirectory$Arbitrary <- generate_unique_strings(length(elementsofinterest))
  tablelookupdirectory$first_write <- TRUE
  
  # Enable progress handler
  #handlers("txtprogressbar") 
  handlers(list(
    handler_progress(
      format   = "lmerEffectsBootstrapSimulationConsolidator() [:bar] :percent :eta",
      width    = 120,
      complete = "="
    )
  ))
  
  file_list <- list.files(results$futuretag$tmpdir, pattern = "^result_.*\\.RData$")
  file_listL <- length(file_list)
  if (file_listL > results$futuretag$repetitions) {
    file_listL <- results$futuretag$repetitions
  }
  
  # load data from files
  resstoreplaceholder <- with_progress({
    p <- progressor(along = 1:file_listL)
    for (i in 1:file_listL) {
      smp <- tryCatch({
        load(file.path(results$futuretag$tmpdir, file_list[i]))  # loads `result`
        # loop through output
        for (cE in 1:nrow(tablelookupdirectory)) {
          textcall <- sprintf('datain <- result$%s', tablelookupdirectory$Real[cE])
          eval(parse(text=textcall))
          
          if (!is.data.frame(datain)) {
            datain <- data.frame(datain)
            if (tablelookupdirectory$Real[cE] == 'numparticipants') {
              colnames(datain) <- c('numparticipants')
            }
          }
          if (tablelookupdirectory$first_write[cE]) {
            save(datain, file = file.path(tmpdir, sprintf(tablelookupdirectory$Arbitrary[cE])))
            tablelookupdirectory$first_write[cE] <- FALSE
            rm(datain)
          } else {
            datainnew <- datain
            load(file = file.path(tmpdir, sprintf(tablelookupdirectory$Arbitrary[cE])))
            datain <- rbind(datain, datainnew)
            save(datain, file = file.path(tmpdir, sprintf(tablelookupdirectory$Arbitrary[cE])))
            rm(datain, datainnew)
          }
        }
        rm(result); gc()
        smp <- TRUE
      }, error = function(e) {
        cat(sprintf('lmerEffectsBootstrap - summary failure\n'))
        smp <- FALSE
      })
      p()
    }
  })
  
  Sys.sleep(1) # to make sure files are read in
  
  # summarize
  tablelookupdirectory$Summarized <- FALSE
  for (cE in 1:nrow(tablelookupdirectory)) {
    load(file = file.path(tmpdir, sprintf(tablelookupdirectory$Arbitrary[cE]))) # load datain
    
    if (tablelookupdirectory$Real[cE] == 'numparticipants') {
      results$meanofparticipants <- mean(unlist(datain$numparticipants, recursive=TRUE), na.rm=TRUE)
      results$medianofparticipants <- median(unlist(datain$numparticipants, recursive=TRUE), na.rm=TRUE)
      results$minofparticipants <- min(unlist(datain$numparticipants, recursive=TRUE), na.rm=TRUE)
      results$maxofparticipants <- max(unlist(datain$numparticipants, recursive=TRUE), na.rm=TRUE)
      tablelookupdirectory$Summarized[cE] <- TRUE
      
    } else if (tablelookupdirectory$Real[cE] == 'descriptives') {
      # descriptives
      textcall <- sprintf('newstats <- results$%s', tablelookupdirectory$Real[cE])
      eval(parse(text=textcall))
      
      for (cR in 1:nrow(newstats)) {
        tempdbs <- datain[which(datain$Group == newstats$Group[cR]),]
        colsofinterest <- c("N", "Missing", "Mean", "Median", "SD", "SE")
        for (cC in 1:length(colsofinterest)) {
          tempvect <- as.numeric(unlist(tempdbs[,colsofinterest[cC]]))
          tempvect <- tempvect[which(!is.na(tempvect))]
          if (average == 'median') {
            newstats[cR,colsofinterest[cC]] <- median(tempvect, na.rm=TRUE)
          } else {
            newstats[cR,colsofinterest[cC]] <- mean(tempvect, na.rm=TRUE)
          }
        }
        distributiontable <- data.frame(matrix(NA, nrow=nrow(tempdbs), ncol=64))
        for (csR in 1:nrow(tempdbs)) {
          textcall <- sprintf('distributiontable[csR,] <- c(%s)', tempdbs$DistributionData[csR])
          eval(parse(text=textcall))
        }
        for (cC in 1:ncol(distributiontable)) {
          distributiontable[,cC] <- as.numeric(distributiontable[,cC])
        }
        if (average == 'median') {
          newstats$DistributionData[cR] <- paste(miscTools::colMedians(distributiontable, na.rm=TRUE), collapse = ",")
        } else {
          newstats$DistributionData[cR] <- paste(colMeans(distributiontable, na.rm=TRUE), collapse = ",")
        }
        decisionindx <- which(tempdbs$DistributionDecision == "Normal")
        if (length(decisionindx) > 0) {
          if ((length(decisionindx) / nrow(tempdbs)) >= 0.6) {
            newstats$DistributionDecision[cR] <- "Normal"
          } else {
            newstats$DistributionDecision[cR] <- "Not Normal"
          }
        }
      }
      
      # dirty but effective
      newstats$Mean <- round(round(round(round(as.numeric(newstats$Mean), digits=4), digits=3), digits=2), digits=1)
      newstats$Median <- round(round(round(round(as.numeric(newstats$Median), digits=4), digits=3), digits=2), digits=1)
      newstats$SD <- round(round(round(round(as.numeric(newstats$SD), digits=4), digits=3), digits=2), digits=1)
      newstats$SE <- round(round(round(round(as.numeric(newstats$SE), digits=4), digits=3), digits=2), digits=1)
      
      newstats$Mean <- sprintf('%.1f', newstats$Mean)
      newstats$Median <- sprintf('%.1f', newstats$Median)
      newstats$SD <- sprintf('%.1f', newstats$SD)
      newstats$SE <- sprintf('%.1f', newstats$SE)
      
      results$descriptives <- newstats
      tablelookupdirectory$Summarized[cE] <- TRUE
      rm(newstats)
      
    } else if ((all(c('SSn', 'SSd', 'F.value') %in% colnames(datain))) | (all(c('LogLikelihood', 'LRT') %in% colnames(datain))) | (all(c('portion', 'ci.lower', 'ci.upper') %in% colnames(datain)))) {
      # fixed, random, r2
      textcall <- sprintf('newstats <- results$%s', tablelookupdirectory$Real[cE])
      eval(parse(text=textcall))
      
      # clear to avoid confusion
      if ('textoutput' %in% colnames(newstats)) {
        newstats$textoutput <- NA
      }
      if ('significance' %in% colnames(newstats)) {
        newstats$significance <- NA
      }
      if ('F.value' %in% colnames(newstats)) {
        newstats$F.value.ci.lower <- NA
        newstats$F.value.ci.upper <- NA
      }
      if ('p.value' %in% colnames(newstats)) {
        newstats$p.value.ci.lower <- NA
        newstats$p.value.ci.upper <- NA
      }
      
      for (cR in 1:nrow(newstats)) {
        if ((all(c('SSn', 'SSd', 'F.value') %in% colnames(datain))) | (all(c('LogLikelihood', 'LRT') %in% colnames(datain)))) {
          tempdbs <- datain[which(datain$Effect == newstats$Effect[cR]),]
        } else if (all(c('portion', 'ci.lower', 'ci.upper') %in% colnames(datain))) {
          tempdbs <- datain[which(datain$portion == newstats$portion[cR]),]
        }
        colsofinterest <- c("DFn", "DFd", "SSn", "SSd", "SSe", "F.value", "DF", "LogLikelihood", "LRT", "p.value", "partialetasquared", "fsquared","fsquared.ci.lower", "fsquared.ci.upper", "effects")
        for (cC in 1:length(colsofinterest)) {
          if (colsofinterest[cC] %in% colnames(tempdbs)) {
            tempvect <- as.numeric(unlist(tempdbs[,colsofinterest[cC]]))
            tempvect <- tempvect[which(!is.na(tempvect))]
            if (average == 'median') {
              newstats[cR,colsofinterest[cC]] <- median(tempvect, na.rm=TRUE)
            } else {
              newstats[cR,colsofinterest[cC]] <- mean(tempvect, na.rm=TRUE)
            }
          }
        }
        
        # significance test
        if (all(c('p.value','significance') %in% colnames(newstats))) {
          outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
          newstats$significance[cR] <- outPvalue$significance
        }
        
        # additional confidence intervals
        if (all(c('F.value') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'F.value']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$F.value[cR]) {ci_boot[1] < newstats$F.value[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$F.value[cR]) {ci_boot[2] < newstats$F.value[cR]} 
          }
          newstats$F.value.ci.lower[cR] <- ci_boot[1]
          newstats$F.value.ci.upper[cR] <- ci_boot[2]
        }
        if (all(c('p.value') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'p.value']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$p.value[cR]) {ci_boot[1] < newstats$p.value[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$p.value[cR]) {ci_boot[2] < newstats$p.value[cR]} 
          }
          newstats$p.value.ci.lower[cR] <- ci_boot[1]
          newstats$p.value.ci.upper[cR] <- ci_boot[2]
        }
        if (all(c('ci.lower', 'ci.upper') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'effects']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$effects[cR]) {ci_boot[1] < newstats$effects[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$effects[cR]) {ci_boot[2] < newstats$effects[cR]} 
          }
          newstats$ci.lower[cR] <- ci_boot[1]
          newstats$ci.upper[cR] <- ci_boot[2]
        }
      }
      # put data back
      textcall <- sprintf('results$%s <- newstats', tablelookupdirectory$Real[cE])
      eval(parse(text=textcall))
      rm(newstats)
      tablelookupdirectory$Summarized[cE] <- TRUE
      
    } else if (all(c('contrast','df', 't.ratio', 'p.value', 'correlation') %in% colnames(datain))) {
      # posthoc
      textcall <- sprintf('newstats <- results$%s', tablelookupdirectory$Real[cE])
      eval(parse(text=textcall))
      
      # clear to avoid confusion
      if ('textoutput' %in% colnames(newstats)) {
        newstats$textoutput <- NA
      }
      if ('significant' %in% colnames(newstats)) {
        newstats$significant <- NA
      }
      if ('t.ratio' %in% colnames(newstats)) {
        newstats$t.conf.int.lower <- NA
        newstats$t.conf.int.upper <- NA
      }
      if ('p.value' %in% colnames(newstats)) {
        newstats$p.conf.int.lower <- NA
        newstats$p.conf.int.upper <- NA
      }
      
      for (cR in 1:nrow(newstats)) {
        tempdbs <- datain[which(datain$idtag == newstats$idtag[cR] & datain$hold == newstats$hold[cR]),]
        if (nrow(tempdbs) == 0) {
          # data simulation likely swapped position of c1 and c2
          tempvect <- stringr::str_split(newstats$idtag[cR], "-")[[1]]
          newtempvect <- tempvect
          newtempvect[length(tempvect)] <- tempvect[length(tempvect)-2]
          newtempvect[length(tempvect)-2] <- tempvect[length(tempvect)]
          newstats$idtag[cR] <- paste0(newtempvect, collapse="-")
          c1name <- newstats$C1name
          c2name <- newstats$C2name
          newstats$C1name <- c2name
          newstats$C2name <- c1name
          
          c1n <- newstats$C1n
          c2n <- newstats$C2n
          newstats$C1n <- c2n
          newstats$C2n <- c1n
          
          c1n <- newstats$C1mean
          c2n <- newstats$C2mean
          newstats$C1mean <- c2n
          newstats$C2mean <- c1n
          
          c1n <- newstats$C1sd
          c2n <- newstats$C2sd
          newstats$C1sd <- c2n
          newstats$C2sd <- c1n
          
          tempdbs <- datain[which(datain$idtag == newstats$idtag[cR] & datain$hold == newstats$hold[cR]),]
        }
        colsofinterest <- c("df", "t.ratio", "p.value", "effectsize", "correlation")
        for (cC in 1:length(colsofinterest)) {
          if (colsofinterest[cC] %in% colnames(tempdbs)) {
            tempvect <- as.numeric(unlist(tempdbs[,colsofinterest[cC]]))
            tempvect <- tempvect[which(!is.na(tempvect))]
            if (average == 'median') {
              newstats[cR,colsofinterest[cC]] <- median(tempvect, na.rm=TRUE)
            } else {
              newstats[cR,colsofinterest[cC]] <- mean(tempvect, na.rm=TRUE)
            }
          }
        }
        
        # significance test
        if (all(c('p.value','significant') %in% colnames(newstats))) {
          outPvalue <- fuzzyP(newstats$p.value[cR], studywiseAlpha=results$studywiseAlpha, html=TRUE)
          newstats$significant[cR] <- outPvalue$significance
        }
        
        # additional confidence intervals
        if (all(c('effectsize') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'effectsize']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$effectsize[cR]) {ci_boot[1] < newstats$effectsize[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$effectsize[cR]) {ci_boot[2] < newstats$effectsize[cR]} 
          }
          newstats$effectsize.conf.int.lower[cR] <- ci_boot[1]
          newstats$effectsize.conf.int.upper[cR] <- ci_boot[2]
        }
        if (all(c('t.ratio') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'t.ratio']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$t.ratio[cR]) {ci_boot[1] < newstats$t.ratio[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$t.ratio[cR]) {ci_boot[2] < newstats$t.ratio[cR]} 
          }
          newstats$t.conf.int.lower[cR] <- ci_boot[1]
          newstats$t.conf.int.upper[cR] <- ci_boot[2]
        }
        if (all(c('p.value') %in% colnames(newstats))) {
          tempvect <- as.numeric(unlist(tempdbs[,'p.value']))
          tempvect <- tempvect[which(!is.na(tempvect))]
          ci_boot <- quantile(tempvect, probs = c(0.05, 0.95), na.rm=TRUE) # one sided
          if (!is.na(ci_boot[1])) {
            if (ci_boot[1] < 0) {ci_boot[1] <- 0}
            if (ci_boot[1] > newstats$p.value[cR]) {ci_boot[1] < newstats$p.value[cR]} 
          }
          if (!is.na(ci_boot[1])) {
            if (ci_boot[2] < 0) {ci_boot[2] <- 0}
            if (ci_boot[2] < newstats$p.value[cR]) {ci_boot[2] < newstats$p.value[cR]} 
          }
          newstats$p.conf.int.lower[cR] <- ci_boot[1]
          newstats$p.conf.int.upper[cR] <- ci_boot[2]
        }
        
        if (reporteddata != 'actual') {
          colsofinterest <- c("C1n", "C1mean", "C1sd", "C2n", "C2mean", "C2sd")
          for (cC in 1:length(colsofinterest)) {
            tempvect <- as.numeric(unlist(tempdbs[,colsofinterest[cC]]))
            tempvect <- tempvect[which(!is.na(tempvect))]
            if (average == 'median') {
              newstats[cR,colsofinterest[cC]] <- median(tempvect, na.rm=TRUE)
            } else {
              newstats[cR,colsofinterest[cC]] <- mean(tempvect, na.rm=TRUE)
            }
          }  
        }
        
      }
      # put data back
      textcall <- sprintf('results$%s <- newstats', tablelookupdirectory$Real[cE])
      eval(parse(text=textcall))
      rm(newstats)
      tablelookupdirectory$Summarized[cE] <- TRUE
    } 
    rm(datain); gc()
  }
  
  # remove folder
  unlink(tmpdir, recursive = TRUE)

  # obtain text outputs
  if ((reporteddata == 'actual') | (reporteddata == 'raw')) {
    subtag <- 'raw'
  } else {
    if (results$futuretag$method == "resample") {
      subtag <- 'resampled'
    } else {
      subtag <- 'simulated'
    }
  }
  results <- lmerEffects2text(results, subtag=subtag, testconfidence=TRUE, significanceconfidence=TRUE)
  
  # see if posthoc adjustments are needed
  if (results$futuretag$boolposthoc) {
    if ((tolower(results$futuretag$methodofposthoccorrection) != 'false') | (tolower(results$futuretag$methodofposthoccorrection) != 'none')) {
      results <- lmerPosthocCorrection(results, method=results$futuretag$methodofposthoccorrection, studywiseAlpha=results$studywiseAlpha, FDRC=0.05)
    }
  }
  
  # tag messageout 
  if ('messageout' %in% names(results)) {
    
    outstring <- sprintf('Unstandardized effects were computed based upon bootstrapped analyses with')
    if (results$futuretag$method == "resample") {
      outstring <- sprintf('%s %d resamples', outstring, file_listL)
      if (!((results$futuretag$resample_min == nrow(stats::model.frame(results$fit))) & (results$futuretag$resample_max == nrow(stats::model.frame(results$fit))))) {
        outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%)',
                             outstring, results$futuretag$resample_min, round((results$futuretag$resample_min/results$futuretag$totalsample)*100, digits=1),
                             resample_max, round((results$futuretag$resample_max/results$futuretag$totalsample)*100, digits=1))
      }
      outstring <- sprintf('%s of the original dataset (with replacement).', outstring)
    }
    if ((results$futuretag$method == "parametric") | (results$futuretag$method == "nonparametric")) {
      if (results$futuretag$method == "parametric") {
        outstring <- sprintf('%s %d datasets simulated from a multivariate normal distribution using the MASS mvrnorm function (Venables &#38; Ripley, 2002).', outstring, file_listL)
      } else {
        outstring <- sprintf('%s %d datasets simulated from a multivariate non-normal distribution using the mnonr unonr function (Qu &#38; Zhang, 2020).', outstring, file_listL)
      }
      outstring <- sprintf('%s For each simulation the covariance matrix was informed by', outstring)
      if (results$futuretag$subsample < 1.0) {
        outstring <- sprintf('%s a subsample of %.1f%% of the original data', outstring, round((results$futuretag$subsample)*100, digits=1))
      } else {
        outstring <- sprintf('%s the full sample of original data', outstring)
      }
      if (results$futuretag$inflation > 1.0) {
        outstring <- sprintf('%s and extrapolated to a final sample of %.0f participants', outstring, round(results$meanofparticipants, digits=0))
      } else {
        outstring <- sprintf('%s.', outstring)
      }
    }
    if (results$futuretag$method == "default") {
      outstring <- sprintf('%s %d datasets simulated by drawing random samples from the conditional distribution of the outcome variable given the estimated model parameters', outstring, file_listL)
      
      if (!((results$futuretag$resample_min == nrow(stats::model.frame(results$fit))) & (results$futuretag$resample_max == nrow(stats::model.frame(results$fit))))) {
        outstring <- sprintf('%s allowing the total number of cases to vary from %d (%.1f%%) to %d (%.1f%%) of the original dataset (with replacement)',
                             outstring, results$futuretag$resample_min, round((results$futuretag$resample_min/results$futuretag$totalsample)*100, digits=1),
                             results$futuretag$resample_max, round((results$futuretag$resample_max/results$futuretag$totalsample)*100, digits=1))
      }
      outstring <- sprintf('%s.', outstring)
      
      if (!is.null(results$futuretag$subsample)) {
        if (results$futuretag$subsample < 1.0) {
          outstring <- sprintf('%s For each simulation the the conditional distribution of the outcome variable was informed by', outstring)
          outstring <- sprintf('%s a subsample of %.1f%% of the original data within each between subjects factor.', outstring, round((results$futuretag$subsample)*100, digits=1))
        }
      }
      
    }
    
    results$messageout <- sprintf('%s %s', results$messageout, outstring)
  }
  
  # remove unnecessary information
  results$futuretag <- NULL
  
  return(results)
}


