#' RmimicTtest
#'
#' @description Compute SPSS style t-test with effect size and confidence intervals. Optional parameters are also provided to compute non-parametric t-tests with appropriate non-parametric effect size estimates. For parametric test the function automatically determines if the variances are equal using levene's test and outputs the correct statistcs.
#'
#' @param data Data frame in long format containing the dependent variable in one column and subject id, between, and within subject variables in other columns.
#' @param dependentvariable Variable name or list of variables to use as dependent variables.
#' @param subjectid Variable name corresponding to the subject id column in the data. Required for paired tests of within subjects variables, encouraged otherwise.
#' @param between Variable name or list of variables to use to compute independent samples t tests.
#' @param within Variable name or list of variables to use to compute paired samples t tests.
#' @param nonparametric Parameter to determine if non parametric tests should be run. Default is FALSE.
#' @param collapse Paramater to determine if multiple observations should be collapsed to a single observation for each participant. Default is TRUE.
#' @param posthoc Parameter to determine what post-hoc comparison approach to use. Options are Bonferroni, Sidak, or Holm-Bonferroni.
#' @param criticaldiff Parameter to specify the critical difference used in Tukey and Scheffe post-hoc comparison approaches.
#' @param confidenceinterval Parameter to control the confidence interval. Default is 0.95.
#' @param studywiseAlpha Parameter to control the study wise alpha for use in the post hoc compairsons. Default is 0.05.
#' @param verbose Parameter to print all output to console. Default is TRUE.
#'
#' @return
#' \item{descriptives}{Data table of descriptive statistics.}
#' \item{stats}{Data table of t test statistics.}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 13, 2019
#'
#' @importFrom stats t.test complete.cases qnorm wilcox.test
#' @importFrom lawstat levene.test
#' @importFrom MBESS conf.limits.nct
#' @importFrom utils packageDate
#' @importFrom pkgcond suppress_conditions
#' 
#' @examples
#' 
#'   ttestresult <- RmimicTtest(PlantGrowth, 
#'           dependentvariable='weight', 
#'           subjectid=NULL, 
#'           between='group', 
#'           within=NULL, 
#'           nonparametric=FALSE,
#'           posthoc="Holm-Bonferroni")
#'
#'
#' @export 

RmimicTtest <- function(data, dependentvariable=NULL, subjectid=NULL, between=NULL, within=NULL, nonparametric=FALSE, collapse=TRUE, posthoc=NULL, criticaldiff=NULL, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE) {
  
  if (!is.null(posthoc)) {
    if (toupper(posthoc) == toupper("Bonferroni")) {
      posthoc <- "Bonferroni"
    } else if (toupper(posthoc) == toupper("Sidak")) {
      posthoc <- "Sidak"
    } else if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
      posthoc <- "Holm-Bonferroni"
    } else if (toupper(posthoc) == toupper("Tukey")) {
      posthoc <- "Tukey"
    } else if (toupper(posthoc) == toupper("Scheffe")) {
      posthoc <- "Scheffe"
    } else if (posthoc == FALSE) {
      posthoc <- NULL
    } else if (toupper(posthoc) == toupper("False Discovery Rate Control")) {
      posthoc <- NULL # handled externally
    } else {
      posthoc <- NULL
    }
  } 
  
  # Assess what got fed into the function
  dependentvariableL <- length(dependentvariable)
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
      stop("Error in RmimicTtest: A paired samples ttest has been requested, but no subjectid parameter was provided to ensure pairwise comparisons. Please provide a column of subject ids.")
    }
  }
  
  # establish output
  masterdescriptives <- NULL
  dataframeout <- data.frame(matrix(NA,nrow=1,ncol=39))
  colnames(dataframeout) <- c('Variable','Comparison','Test','statistic', 'df', 'p.value', 'conf.int.lower', 'conf.int.upper', 'alternative', 'method', 'z.value', 'effectsize', 'effectsize.conf.int.lower', 'effectsize.conf.int.upper', 'correlation', 'correlation.p.value', 'stud.conf.int', 'textoutput', 'Group1.label', 'Group1.n', 'Group1.missing', 'Group1.mean', 'Group1.median', 'Group1.sd', 'Group1.se', 'Group1.min', 'Group1.max', 'Group1.distribution', 'Group2.label', 'Group2.n', 'Group2.missing', 'Group2.mean', 'Group2.median', 'Group2.sd', 'Group2.se', 'Group2.min', 'Group2.max', 'Group2.distribution', 'TestChunk')
  dataframeoutL <- 1
  chunkingL <- 1
  
  # Address Independent Samples tests first
  if (betweenvariableL > 0) {
    # for each between subjects variable
    for (cB in 1:betweenvariableL) {
      
      # subset data for only that variable
      workingdatabase <- data[,c(subjectid[1], between[cB], dependentvariable)]
      
      # populate list of comparisons
      spfactors <- sort(unique(unlist(as.character(data[,between[cB]]))))
      outlist <- Rmimic::determineallpossiblecombinations(spfactors)
      spfactorscomparisons <- c()
      for (cC in 1:length(outlist)) {
        if (lengths(regmatches(outlist[cC], gregexpr(":", outlist[cC]))) == 1) {
          spfactorscomparisons <- c(spfactorscomparisons, outlist[cC])
        }
      }
      rm(cC, outlist, spfactors)
      
      # for each dependent variable
      for (cDV in 1:dependentvariableL) {
        
        # for each comparison
        for (cComparison in 1:length(spfactorscomparisons)) {
          tempvect <- unlist(strsplit(spfactorscomparisons[cComparison], ':'))
          # Obtain databases for each group
          group1data <- workingdatabase[which(workingdatabase[,between[cB]] == tempvect[1]),]
          group2data <- workingdatabase[which(workingdatabase[,between[cB]] == tempvect[2]),]
            
          # Collapse to a single data point for each participant (if requested)
          if (collapse != FALSE) {
            newgroup1data <- data.frame(matrix(NA,nrow=1,ncol=ncol(group1data)))
            colnames(newgroup1data) <- colnames(group1data)
            ids <- unlist(as.character(group1data[,subjectid[1]]))
            for (cids in 1:length(ids)) {
              # subset only data for that participant
              subworkingdatabase <- group1data[which(group1data[,subjectid[1]] == ids[cids]),]
              newgroup1data[cids,] <- subworkingdatabase[1,]
              newgroup1data[cids,dependentvariable[cDV]] <- mean(subworkingdatabase[,dependentvariable[cDV]], na.rm= TRUE)
            }
            group1data <- newgroup1data
            
            newgroup2data <- data.frame(matrix(NA,nrow=1,ncol=ncol(group2data)))
            colnames(newgroup2data) <- colnames(group2data)
            ids <- unlist(as.character(group2data[,subjectid[1]]))
            for (cids in 1:length(ids)) {
              # subset only data for that participant
              subworkingdatabase <- group2data[which(group2data[,subjectid[1]] == ids[cids]),]
              newgroup2data[cids,] <- subworkingdatabase[1,]
              newgroup2data[cids,dependentvariable[cDV]] <- mean(subworkingdatabase[,dependentvariable[cDV]], na.rm= TRUE)
            }
            group2data <- newgroup2data
            rm(newgroup1data, newgroup2data, cids, subworkingdatabase)
          }
        
          # populate data
          comparison1 <- group1data[,dependentvariable[cDV]]
          comparison2 <- group2data[,dependentvariable[cDV]]
          tempframe <- data.frame(DV = c(comparison1, comparison2), Group = c(rep_len("1",length(comparison1)), rep_len("2",length(comparison2))))
          
          # obtain descriptives
          desc <- Rmimic::descriptives(variables='DV', groupvariable='Group', data=tempframe, verbose=FALSE)
          tempframe <- tempframe[stats::complete.cases(tempframe),]
          comparison1 <- comparison1[stats::complete.cases(comparison1)]
          comparison2 <- comparison2[stats::complete.cases(comparison2)]
          
          # determine type of independent test to do
          if (nonparametric == FALSE) {
            # independent samples parametric test
            
            # test variance
            varianceEqual <- TRUE
            variancetest <- pkgcond::suppress_conditions(lawstat::levene.test(tempframe[,'DV'], tempframe[,'Group'], location="median"))
            if (variancetest$p.value <= studywiseAlpha) {
              varianceEqual <- FALSE
            }
            ttestresult <- stats::t.test(x=comparison1, y=comparison2, alternative='two.sided', paired=FALSE, var.equal=varianceEqual, conf.level=confidenceinterval)
            ttestresult$statistic <- abs(ttestresult$statistic)
            ttestresult$effectsize <- abs(ttestresult$statistic * sqrt((1/length(comparison1)) + (1/length(comparison2))))
            temptstat <- ttestresult$statistic
            if (temptstat > 37.6) {
              temptstat <- temptstat
            }
            ncp <- pkgcond::suppress_conditions(MBESS::conf.limits.nct(ncp = temptstat, df = ttestresult$parameter, conf.level = confidenceinterval))
            ttestresult$effectsize.conf.int.lower <- ncp$Lower.Limit * sqrt((1/length(comparison1)) + (1/length(comparison2)))
            ttestresult$effectsize.conf.int.upper <- ncp$Upper.Limit * sqrt((1/length(comparison1)) + (1/length(comparison2)))
            ttestresult$stud.conf.int <- confidenceinterval
            tempout <- Rmimic::ttest2text(ttestresult, verbose=FALSE)
            
            dataframeout[dataframeoutL,'Test'] <- 'Independent samples t-test'
            rm(varianceEqual, variancetest, ncp, ttestresult)
            
          } else {
            # mann-whitney t-test (independent samples)
            
            # The U statistic changes depending on the order the variables are entered,
            # Convention says take the smallest U value as the P value remains the same either way.
            ttest1 <- stats::wilcox.test(x=comparison1, y=comparison2, paired=FALSE, correct=FALSE, exact=FALSE, conf.int=TRUE, conf.level=confidenceinterval)
            ttest2 <- stats::wilcox.test(x=comparison2, y=comparison1, paired=FALSE, correct=FALSE, exact=FALSE, conf.int=TRUE, conf.level=confidenceinterval)
            if (ttest1$statistic < ttest2$statistic) {
              ttestresult <- ttest1
              ttestRev <- ttest2
            } else {
              ttestresult <- ttest2
              ttestRev <- ttest1
            }
            ttestresult$z.value <- abs(stats::qnorm(ttestresult$p.value/ 2))
            ttestresult$effectsize <- abs(stats::qnorm(ttestresult$p.value/ 2))/sqrt(length(comparison1)+length(comparison2))
            meandiff <- mean(comparison2)-mean(comparison1)
            if (!(ttestresult$conf.int[1] <= meandiff) & (ttestresult$conf.int[2] >= meandiff)) {
              ttestresult$conf.int[1] <- ttestRev$conf.int[1]
              ttestresult$conf.int[2] <- ttestRev$conf.int[2]
            }
            tempout <- Rmimic::ttest2text(ttestresult, verbose=FALSE)
            
            dataframeout[dataframeoutL,'Test'] <- 'Independent samples Mann-Whitney test'
            rm(ttest1, ttest2, ttestRev, meandiff, ttestresult)
          }
          
          # populate output
          dataframeout[dataframeoutL,'Variable'] <- dependentvariable[cDV]
          dataframeout[dataframeoutL,'Comparison'] <- sprintf('%s: %s-%s', between[cB], tempvect[1], tempvect[2])
          dataframeout[dataframeoutL,names(tempout$statistics)] <- tempout$statistics
          dataframeout[dataframeoutL,'textoutput'] <- tempout$text
          dataframeout[dataframeoutL,'Group1.label'] <- tempvect[1]
          dataframeout[dataframeoutL,c('Group1.n', 'Group1.missing', 'Group1.mean', 'Group1.median', 'Group1.sd', 'Group1.se', 'Group1.min', 'Group1.max', 'Group1.distribution')] <- desc[1,c('N', 'Missing', 'Mean', "Median", 'SD', 'SE', 'Min', 'Max', 'Distribution')]
          dataframeout[dataframeoutL,'Group2.label'] <- tempvect[2]
          dataframeout[dataframeoutL,c('Group2.n', 'Group2.missing', 'Group2.mean', 'Group2.median', 'Group2.sd', 'Group2.se', 'Group2.min', 'Group2.max', 'Group2.distribution')] <- desc[2,c('N', 'Missing', 'Mean', "Median", 'SD', 'SE', 'Min', 'Max', 'Distribution')]
          dataframeout[dataframeoutL,'TestChunk'] <- chunkingL
          
          desc$Variable <- dependentvariable[cDV]
          desc$Group <- between[cB]
          desc$CollapsedName[1] <- sprintf('%s: %s', between[cB], tempvect[1])
          desc$CollapsedName[2] <- sprintf('%s: %s', between[cB], tempvect[2])
          if (!is.null(masterdescriptives)) {
            masterdescriptives <- rbind(masterdescriptives, desc)
          } else {
            masterdescriptives <- desc
          }
          
          dataframeoutL <- dataframeoutL + 1
          rm(tempout, tempframe, desc, comparison1, comparison2)
          
          
        } # end cComparison
        chunkingL <- chunkingL + 1
      } # end cDV
      
      rm(cComparison, workingdatabase, spfactorscomparisons)
    } # end cB
    
    rm(cB)
  } # end betweenvariableL
  
  # Address Paired Samples tests next
  if (withinvariableL > 0) {
    # for each between subjects variable
    for (cB in 1:withinvariableL) {
      
      # subset data for only that variable
      workingdatabase <- data[,c(subjectid[1], within[cB], dependentvariable)]
      
      # populate list of comparisons
      spfactors <- sort(unique(unlist(as.character(data[,within[cB]]))))
      outlist <- Rmimic::determineallpossiblecombinations(spfactors)
      spfactorscomparisons <- c()
      for (cC in 1:length(outlist)) {
        if (lengths(regmatches(outlist[cC], gregexpr(":", outlist[cC]))) == 1) {
          spfactorscomparisons <- c(spfactorscomparisons, outlist[cC])
        }
      }
      
      # for each dependent variable
      for (cDV in 1:dependentvariableL) {
      
        # for each comparison
        for (cComparison in 1:length(spfactorscomparisons)) {
          tempvect <- unlist(strsplit(spfactorscomparisons[cComparison], ':'))
          # Obtain databases for each group
          group1data <- workingdatabase[which(workingdatabase[,within[cB]] == tempvect[1]),]
          group2data <- workingdatabase[which(workingdatabase[,within[cB]] == tempvect[2]),]
          
          # Make sure each participants data has a match - paired samples
          subgroup1data <- group1data[,c(subjectid[1], dependentvariable[cDV])]
          subgroup1data <- subgroup1data[stats::complete.cases(subgroup1data),]
          subgroup2data <- group2data[,c(subjectid[1], dependentvariable[cDV])]
          subgroup2data <- subgroup2data[stats::complete.cases(subgroup2data),]
          
          if (collapse != FALSE) {
            # Collapse to a single data point for each participant (if requested)
            
            newsubgroup1data <- data.frame(matrix(NA,nrow=1,ncol=ncol(subgroup1data)))
            colnames(newsubgroup1data) <- colnames(subgroup1data)
            ids <- unlist(as.character(subgroup1data[,subjectid[1]]))
            for (cids in 1:length(ids)) {
              # subset only data for that participant
              subworkingdatabase <- subgroup1data[which(subgroup1data[,subjectid[1]] == ids[cids]),]
              newsubgroup1data[cids,] <- subworkingdatabase[1,]
              newsubgroup1data[cids,dependentvariable[cDV]] <- mean(subworkingdatabase[,dependentvariable[cDV]], na.rm= TRUE)
            }
            subgroup1data <- newsubgroup1data
            
            newsubgroup2data <- data.frame(matrix(NA,nrow=1,ncol=ncol(subgroup2data)))
            colnames(newsubgroup2data) <- colnames(subgroup2data)
            ids <- unlist(as.character(subgroup2data[,subjectid[1]]))
            for (cids in 1:length(ids)) {
              # subset only data for that participant
              subworkingdatabase <- subgroup2data[which(subgroup2data[,subjectid[1]] == ids[cids]),]
              newsubgroup2data[cids,] <- subworkingdatabase[1,]
              newsubgroup2data[cids,dependentvariable[cDV]] <- mean(subworkingdatabase[,dependentvariable[cDV]], na.rm= TRUE)
            }
            subgroup2data <- newsubgroup2data
            rm(newsubgroup1data, newsubgroup2data, cids, subworkingdatabase)
          }
          tempcompdatabase <- merge(subgroup1data, subgroup2data, by=subjectid[1], all = FALSE)
          
          # populate data
          comparison1 <- tempcompdatabase[,2]
          comparison2 <- tempcompdatabase[,3]
          rm(tempcompdatabase, subgroup1data, subgroup2data)
          tempframe <- data.frame(DV = c(comparison1, comparison2), Group = c(rep_len("1",length(comparison1)), rep_len("2",length(comparison2))))
          
          # obtain descriptives
          desc <- Rmimic::descriptives(variables='DV', groupvariable='Group', data=tempframe, verbose=FALSE)
          
          # determine type of paired test to do
          if (nonparametric == FALSE) {
            # paired samples parametric test
            
            # test variance
            varianceEqual <- TRUE
            variancetest <- lawstat::levene.test(tempframe[,'DV'], tempframe[,'Group'], location="median")
            if (variancetest$p.value <= studywiseAlpha) {
              varianceEqual <- FALSE
            }
            ttestresult <- stats::t.test(x=comparison1, y=comparison2, alternative='two.sided', paired=TRUE, var.equal=varianceEqual, conf.level=confidenceinterval)
            ttestresult$statistic <- abs(ttestresult$statistic)
            correlationtest <- stats::cor.test(comparison1, comparison2, alternative='two.sided', method = "pearson", conf.level = confidenceinterval, use = "complete.obs")
            ttestresult$correlation <- correlationtest$estimate[[1]]
            ttestresult$correlation.p.value <- correlationtest$p.value[[1]]
            ttestresult$effectsize <- ttestresult$statistic * sqrt((2*(1-correlationtest$estimate[[1]]))/length(comparison2))
            temptstat <- ttestresult$statistic
            if (temptstat > 37.6) {
              temptstat <- temptstat
            }
            ncp <- pkgcond::suppress_conditions(MBESS::conf.limits.nct(ncp = temptstat, df = ttestresult$parameter, conf.level = confidenceinterval))
            ttestresult$effectsize.conf.int.lower <- ncp$Lower.Limit * sqrt((2*(1-correlationtest$estimate[[1]]))/length(comparison2))
            ttestresult$effectsize.conf.int.upper <- ncp$Upper.Limit * sqrt((2*(1-correlationtest$estimate[[1]]))/length(comparison2))
            ttestresult$stud.conf.int <- confidenceinterval
            tempout <- Rmimic::ttest2text(ttestresult, verbose=FALSE)
            
            dataframeout[dataframeoutL,'Test'] <- 'Paired samples t-test'
            rm(varianceEqual, variancetest, correlationtest, ncp, ttestresult)
            
          } else {
            #Wilcoxon signed rank test (paired samples)
            
            # The V statistic changes depending on the order the variables are entered,
            # Convention says take the smallest V value as the P value remains the same either way.
            ttest1 <- stats::wilcox.test(x=comparison1, y=comparison2, paired=TRUE, correct=FALSE, exact=FALSE, conf.int=TRUE, conf.level=confidenceinterval)
            ttest2 <- stats::wilcox.test(x=comparison2, y=comparison1, paired=TRUE, correct=FALSE, exact=FALSE, conf.int=TRUE, conf.level=confidenceinterval)
            if (ttest1$statistic < ttest2$statistic) {
              ttestresult <- ttest1
              ttestRev <- ttest2
            } else {
              ttestresult <- ttest2
              ttestRev <- ttest1
            }
            ttestresult$z.value <- abs(stats::qnorm(ttestresult$p.value/ 2))
            ttestresult$effectsize <- abs(stats::qnorm(ttestresult$p.value/ 2))/sqrt(length(comparison2))
            meandiff <- mean(comparison2)-mean(comparison1)
            if (!(ttestresult$conf.int[1] <= meandiff) & (ttestresult$conf.int[2] >= meandiff)) {
              ttestresult$conf.int[1] <- ttestRev$conf.int[1]
              ttestresult$conf.int[2] <- ttestRev$conf.int[2]
            }
            tempout <- Rmimic::ttest2text(ttestresult, verbose=FALSE)
            
            dataframeout[dataframeoutL,'Test'] <- 'Paired samples Wilcoxon signed rank test'
            rm(ttest1, ttest2, ttestRev, meandiff, ttestresult)
            
          }
          
          # populate output
          dataframeout[dataframeoutL,'Variable'] <- dependentvariable[cDV]
          dataframeout[dataframeoutL,'Comparison'] <- sprintf('%s: %s-%s', within[cB], tempvect[1], tempvect[2])
          dataframeout[dataframeoutL,names(tempout$statistics)] <- tempout$statistics
          dataframeout[dataframeoutL,'textoutput'] <- tempout$text
          dataframeout[dataframeoutL,'Group1.label'] <- tempvect[1]
          dataframeout[dataframeoutL,c('Group1.n', 'Group1.missing', 'Group1.mean', 'Group1.median', 'Group1.sd', 'Group1.se', 'Group1.min', 'Group1.max', 'Group1.distribution')] <- desc[1,c('N', 'Missing', 'Mean', "Median", 'SD', 'SE', 'Min', 'Max', 'Distribution')]
          dataframeout[dataframeoutL,'Group2.label'] <- tempvect[2]
          dataframeout[dataframeoutL,c('Group2.n', 'Group2.missing', 'Group2.mean', 'Group2.median', 'Group2.sd', 'Group2.se', 'Group2.min', 'Group2.max', 'Group2.distribution')] <- desc[2,c('N', 'Missing', 'Mean', "Median", 'SD', 'SE', 'Min', 'Max', 'Distribution')]
          dataframeout[dataframeoutL,'TestChunk'] <- chunkingL
          
          desc$Variable <- dependentvariable[cDV]
          desc$Group <- within[cB]
          desc$CollapsedName[1] <- sprintf('%s: %s', within[cB], tempvect[1])
          desc$CollapsedName[2] <- sprintf('%s: %s', within[cB], tempvect[2])
          if (!is.null(masterdescriptives)) {
            masterdescriptives <- base::rbind(masterdescriptives, desc)
          } else {
            masterdescriptives <- desc
          }
          
          dataframeoutL <- dataframeoutL + 1
          rm(tempout, tempframe, desc, comparison1, comparison2)
          
        } # end cComparison
        chunkingL <- chunkingL + 1
      } # end cDV
    } # end cB
  } # end withinvariableL
  
  # remove duplicate demographic entries
  finalmasterdescriptives <- data.frame(matrix(NA, nrow=1, ncol=ncol(masterdescriptives)))
  finalmasterdescriptivesL <- 1
  colnames(finalmasterdescriptives) <- colnames(masterdescriptives)
  mdVar <- unique(masterdescriptives$Variable)
  for (cDV in 1:length(mdVar)) {
    tempdataframe <- masterdescriptives[which(masterdescriptives$Variable == mdVar[cDV]),]
    mdGro <- unique(tempdataframe$Group)
    for (cB in 1:length(mdGro)) {
      subtempdataframe <- tempdataframe[which(tempdataframe$Group == mdGro[cB]),]
      mdNam <- unique(subtempdataframe$CollapsedName)
      for (cN in 1:length(mdNam)) {
        firstindex <- which(subtempdataframe$CollapsedName == mdNam[cN])[1]
        finalmasterdescriptives[finalmasterdescriptivesL,] <- subtempdataframe[firstindex,]
        finalmasterdescriptivesL <- finalmasterdescriptivesL + 1
      }
    }
  }
  rm(masterdescriptives, subtempdataframe, tempdataframe, mdNam, mdGro, mdVar, cN, cB, cDV)
  
  # populate interpretation
  dataframeout$interpretation <- NA
  for (cI in 1:nrow(dataframeout)) {
    temptextout0 <- "No significant differences were observed between"
    temptextout1 <- sprintf("%s (%.1f +/- %.1f) and %s (%.1f +/- %.1f)", dataframeout$Group1.label[cI], dataframeout$Group1.mean[cI], dataframeout$Group1.sd[cI], dataframeout$Group2.label[cI], dataframeout$Group2.mean[cI], dataframeout$Group2.sd[cI])
    temptextout2 <- ""
    if (dependentvariableL > 1) {
      temptextout3 <- sprintf(" for %s; %s", dataframeout$Variable[cI], dataframeout$textoutput[cI])
    } else {
      temptextout3 <- sprintf("; %s", dataframeout$textoutput[cI])
    }
    
    outPvalue <- fuzzyP(dataframeout$p.value[cI])
    if (outPvalue$interpret <= studywiseAlpha) {
      temptextout0 <- "The difference between"
      temptextout2 <- "was statistically significant"
      dataframeout$interpretation[cI] <- sprintf("%s %s %s%s", temptextout0, temptextout1, temptextout2, temptextout3)
    } else {
      dataframeout$interpretation[cI] <- sprintf("%s %s%s", temptextout0, temptextout1, temptextout3)
    }
    rm(outPvalue)
  }
  rm(cI, temptextout0, temptextout1, temptextout2, temptextout3)
  
  # Post hoc adjustments
  if (!is.null(posthoc)) {
    chunklabels <- unique(as.character(dataframeout$TestChunk))
    for (cC in 1:length(chunklabels)) {
      submatrix <- dataframeout[which(dataframeout$TestChunk == chunklabels[cC]),] # subset for test chunk
      
      if ((toupper(posthoc) == toupper("Bonferroni")) | (toupper(posthoc) == toupper("Sidak"))) {
        critp <- studywiseAlpha
        criticalphrase <- ""
        if (toupper(posthoc) == toupper("Bonferroni")) {
          critp <- (studywiseAlpha/nrow(submatrix))
          criticalphrase <- sprintf("(Bonferroni critical alpha = %.4f)", round(critp,4))
        } else if (toupper(posthoc) == toupper("Sidak")) {
          critp <- (1-(1-studywiseAlpha)^(1/nrow(submatrix)))
          criticalphrase <- sprintf("(Sidak critical alpha = %.4f)", round(critp,4))
        }
        for (cR in 1:nrow(submatrix)) {
          outPvalue <- fuzzyP(submatrix$p.value[cR])
          if (outPvalue$interpret <= studywiseAlpha) {
            if (outPvalue$interpret > critp) {
              submatrix$interpretation[cR] <- sprintf('%s However, that difference did not remain significant following correction for multiple comparisons, %s.', submatrix$interpretation[cR], criticalphrase)
            }
          }
        }
      }
      
      if ((toupper(posthoc) == toupper("Tukey")) | (toupper(posthoc) == toupper("Scheffe"))) {
        if (!is.null(criticaldiff)) {
          for (cR in 1:nrow(dataframeout)) {
            actdiff <- abs(as.numeric(dataframeout$Group1.mean[cR]) - as.numeric(dataframeout$Group2.mean[cR]))
            if (actdiff < as.numeric(criticaldiff)) {
              if (toupper(posthoc) == toupper("Tukey")) {
                criticalphrase <- sprintf("(Tukey critical difference = %.2f)", round(criticaldiff,digits=2))
              } else {
                criticalphrase <- sprintf("(Scheffe critical difference = %.2f)", round(criticaldiff,digits=2))
              }
              dataframeout$interpretation[cR] <- sprintf('%s However, that difference did not remain significant following correction for multiple comparisons, %s.', dataframeout$interpretation[cR], criticalphrase)
            }
          }
        }
      }
      
      if (toupper(posthoc) == toupper("Holm-Bonferroni")) {
        #Holm, S. 1979. A simple sequential rejective multiple test procedure. Scandinavian Journal of Statistics 6:65-70
        temp <- data.frame(matrix(NA, nrow=nrow(submatrix), ncol=2))
        colnames(temp) <- c("p.value", "location")
        for (cR in 1:nrow(submatrix)) {
          outPvalue <- fuzzyP(submatrix$p.value[cR])
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
                submatrix$interpretation[temp$location[rank]] <- sprintf('%s However, that difference did not remain significant following correction for multiple comparisons, %s.', submatrix$interpretation[temp$location[rank]], criticalphrase)
              }
            } else {
              # P value is not significant
              rank <- rank + 1
            }
          }
          rank <- rank + 1
        }
      }
      dataframeout[which(dataframeout$TestChunk == chunklabels[cC]),] <- submatrix # return values
    }
  }
  
  # output to console
  if (verbose == TRUE) {
    
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
    
    temptext <- "T-test Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- ""
    if (betweenvariableL > 0) {
      #independent samples ttest
      outstring <- sprintf('%sBetween subject analysis were conducted using', outstring)
      if (nonparametric == FALSE) {
        outstring <- sprintf('%s an independent samples t-test', outstring)
      } else {
        outstring <- sprintf('%s the non-parametric Mann-Whitney U-test', outstring)
      }
      if ((!is.null(posthoc)) & (posthoc != FALSE)) {
        outstring <- sprintf('%s using the %s approach for post-hoc comparison corrections', outstring, posthoc)
      }
      outstring <- sprintf('%s.', outstring)
    }
    if (withinvariableL > 0) {
      outstring <- sprintf('%sWithin subject analysis were conducted using', outstring)
      if (nonparametric == FALSE) {
        outstring <- sprintf('%s a paired samples t-test', outstring)
      } else {
        outstring <- sprintf('%s the non-parametric Wilcoxon signed rank test', outstring)
      }
      if ((!is.null(posthoc)) & (posthoc != FALSE)) {
        outstring <- sprintf('%s using the %s approach for post-hoc comparison corrections', outstring, posthoc)
      }
      outstring <- sprintf('%s.', outstring)
    }
      
    outstring <- sprintf('%s Analysis were conducted using the', outstring)
    outstring <- sprintf('%s stats (R Core Team, %s)', outstring, strsplit(as.character(utils::packageDate("stats")),"-")[[1]][1])
    outstring <- sprintf('%s, MBESS (Kelley, %s)', outstring, strsplit(as.character(utils::packageDate("MBESS")),"-")[[1]][1])
    outstring <- sprintf('%s, and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    
    # find out what needs to be printed
    grouplabels <- unique(as.character(finalmasterdescriptives$CollapsedName))
    nongrouplabels <- unique(as.character(finalmasterdescriptives$Variable))
    ngroups <- length(grouplabels)
    if (ngroups > 1) {
      if (sum(unlist(finalmasterdescriptives$Missing), na.rm = TRUE) > 0) {
        vectnames <- c("Variable", "CollapsedName", "N", "Missing", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      } else {
        vectnames <- c("Variable", "CollapsedName", "N", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      }
    } else {
      if (sum(unlist(finalmasterdescriptives$Missing), na.rm = TRUE) > 0) {
        vectnames <- c("Variable", "N", "Missing", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      } else {
        vectnames <- c("Variable", "N", "Mean", "Median", "SD", "SE", "Min", "Max", "Distribution")
      }
    }
    sepgap <- matrix(floor(spansize/length(vectnames)), nrow=1, ncol=length(vectnames))
    colnames(sepgap) <- vectnames
    sepgap[1,1] <- sepgap[1,1] + 2
    
    # Write header
    summarylabel <- "Descriptive Statistics"
    cat(sprintf("\n%s\n", summarylabel))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    for (cC in 1:length(vectnames)) {
      if (vectnames[cC] == "CollapsedName") {
        cat(sprintf("%-*s",sepgap[1,cC],"Group"))
      } else {
        cat(sprintf("%-*s",sepgap[1,cC],vectnames[cC]))
      }
    }
    cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    for (cB in 1:length(nongrouplabels)) {
      # Write Variable Label
      cat(sprintf("%s\n",nongrouplabels[cB]))
      for (cG in 1:length(grouplabels)) {
        submatrix <- finalmasterdescriptives[which(finalmasterdescriptives$Variable == nongrouplabels[cB]),] # subset for variable
        submatrix <- submatrix[which(submatrix$CollapsedName == grouplabels[cG]),] # subset for group
        cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:1]))), " "), collapse = "")))
        starval <- 2
        if (vectnames[2] == "CollapsedName") {
          # Indent properly
          cat(sprintf("%s\n",grouplabels[cG]))
          cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:2]))), " "), collapse = "")))
          starval <- 3
        }
        
        # Cycle Through Values
        for (cC in starval:length(vectnames)) {
          if ((is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
            # value is numeric and not missing
            pullvalue <- submatrix[1,sprintf("%s", vectnames[cC])]
            if (is.integer(pullvalue)) {
              pullvalue <- sprintf('%d', pullvalue)
            } else {
              pullvalue <- sprintf('%0.1f', round(pullvalue, digits=1))
            }
            cat(sprintf("%-*s",sepgap[1,cC],pullvalue))
          } else {
            if ((!is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
              cat(sprintf("%-*s",sepgap[1,cC],submatrix[1,sprintf("%s", vectnames[cC])]))
            } else {
              cat(sprintf("%-*s",sepgap[1,cC],"na"))
            }
          }
        }
        cat(sprintf("\n"))
      }
      if (cB < length(nongrouplabels)) {
        cat(sprintf("%s\n",paste(replicate(floor(spansize/3)-1, bigspancharacter), collapse = "")))
      }
    }
    cat(sprintf("%s\n\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    testlabels <- unique(as.character(dataframeout$Test))
    testlabelsL <- length(testlabels)
    for (cT in 1:testlabelsL) {
      sepgap <- matrix(floor(spansize/6), nrow=1, ncol=6)
      if (testlabels[cT] == "Independent samples t-test") {
        vectnames <- c("Variable", "CollapsedName", "df", "statistic", "p.value", "effectsize")
        sepgap[1,3] <- sepgap[1,3] - 4
        sepgap[1,4] <- sepgap[1,4] - 4
        sepgap[1,5] <- sepgap[1,5] - 3
      } else if (testlabels[cT] == 'Independent samples Mann-Whitney test') {
        vectnames <- c("Variable", "CollapsedName", "statistic", "z.value", "p.value", "effectsize")
        sepgap[1,2] <- sepgap[1,2] + 6
      } else if (testlabels[cT] == 'Paired samples t-test') {
        vectnames <- c("Variable", "CollapsedName", "df", "statistic", "p.value", "effectsize")
        sepgap[1,3] <- sepgap[1,3] - 4
        sepgap[1,4] <- sepgap[1,4] - 4
        sepgap[1,5] <- sepgap[1,5] - 3
      } else if (testlabels[cT] == 'Paired samples Wilcoxon signed rank test') {
        vectnames <- c("Variable", "CollapsedName", "statistic", "z.value", "p.value", "effectsize")
        sepgap[1,2] <- sepgap[1,2] + 6
      }
      sepgap[1,1] <- sepgap[1,1] - 2
      colnames(sepgap) <- vectnames
      
      # Write header
      cat(sprintf("\n%s\n", testlabels[cT]))
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      for (cC in 1:length(vectnames)) {
        if (vectnames[cC] == "CollapsedName") {
          cat(sprintf("%-*s",sepgap[1,cC],"Group"))
        } else if (vectnames[cC] == "statistic") {
          if (testlabels[cT] == "Independent samples t-test") {
            cat(sprintf("%-*s",sepgap[1,cC],"t"))
          } else if (testlabels[cT] == 'Independent samples Mann-Whitney test') {
            cat(sprintf("%-*s",sepgap[1,cC],"U"))
          } else if (testlabels[cT] == 'Paired samples t-test') {
            cat(sprintf("%-*s",sepgap[1,cC],"t"))
          } else if (testlabels[cT] == 'Paired samples Wilcoxon signed rank test') {
            cat(sprintf("%-*s",sepgap[1,cC],"V"))
          }
        } else if (vectnames[cC] == "z.value") {
          cat(sprintf("%-*s",sepgap[1,cC],"Z"))
        } else if (vectnames[cC] == "p.value") {
          cat(sprintf("%-*s",sepgap[1,cC],"p"))
        } else if (vectnames[cC] == "effectsize") {
          if (testlabels[cT] == "Independent samples t-test") {
            cat(sprintf("%-*s",sepgap[1,cC],sprintf('cohens ds [%d%% CI]', floor(confidenceinterval * 100))))
          } else if (testlabels[cT] == 'Independent samples Mann-Whitney test') {
            cat(sprintf("%-*s",sepgap[1,cC],"r"))
          } else if (testlabels[cT] == 'Paired samples t-test') {
            cat(sprintf("%-*s",sepgap[1,cC],sprintf('cohens drm [%d%% CI]', floor(confidenceinterval * 100))))
          } else if (testlabels[cT] == 'Paired samples Wilcoxon signed rank test') {
            cat(sprintf("%-*s",sepgap[1,cC],"r"))
          }
        } else {
          cat(sprintf("%-*s",sepgap[1,cC],vectnames[cC]))
        }
      }
      cat(sprintf("\n%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
      # output data
      nongrouplabels <- unique(as.character(dataframeout$Variable))
      grouplabels <- unique(as.character(dataframeout$Comparison))
      ngroups <- length(grouplabels)
      for (cB in 1:length(nongrouplabels)) {
        # Write Variable Label
        cat(sprintf("%s\n",nongrouplabels[cB]))
        for (cG in 1:length(grouplabels)) {
          submatrix <- dataframeout[which(dataframeout$Test == testlabels[cT]),] # subset for test
          submatrix <- submatrix[which(submatrix$Variable == nongrouplabels[cB]),] # subset for variable
          submatrix <- submatrix[which(submatrix$Comparison == grouplabels[cG]),] # subset for comparison
          cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:1]))), " "), collapse = "")))
          starval <- 2
          if (vectnames[2] == "CollapsedName") {
            # Indent properly
            cat(sprintf("%s\n",grouplabels[cG]))
            cat(sprintf("%s",paste(replicate(as.integer(sum(unlist(sepgap[1,1:2]))), " "), collapse = "")))
            starval <- 3
          }
          
          # Cycle Through Values
          for (cC in starval:length(vectnames)) {
            if ((is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
              # value is numeric and not missing
              pullvalue <- submatrix[1,sprintf("%s", vectnames[cC])]
              
              # control number of digits
              if ((vectnames[cC] == "statistic") | (vectnames[cC] == "z.value")) {
                pullvalue <- sprintf('%0.1f', round(pullvalue, digits=1))
              } else if (vectnames[cC] == "df") {
                if (is.integer(pullvalue)) {
                  pullvalue <- sprintf('%d', pullvalue)
                } else {
                  pullvalue <- sprintf('%0.1f', round(pullvalue, digits=1))
                }
                if (substr(pullvalue, nchar(pullvalue), nchar(pullvalue)) == "0") {
                  if (substr(pullvalue, nchar(pullvalue)-1, nchar(pullvalue)-1) == ".") {
                    pullvalue <- substr(pullvalue, 1, nchar(pullvalue)-2)
                  }
                }
              } else if (vectnames[cC] == "effectsize") {
                pullvalue <- sprintf('%0.2f', round(pullvalue, digits=2))
                if ((testlabels[cT] == "Independent samples t-test") | (testlabels[cT] == 'Paired samples t-test')) {
                  pullvalue <- sprintf('%s [%d%% CI:', pullvalue, floor(submatrix[[1,'stud.conf.int']] * 100))
                  pullvalue <- sprintf('%s %s to %s]', pullvalue, sprintf('%.2f', round(as.double(submatrix[[1,'effectsize.conf.int.lower']]), digits = 2)), sprintf('%.2f', round(as.double(submatrix[[1,'effectsize.conf.int.upper']]), digits = 2)))
                }
              } else if (vectnames[cC] == "p.value") {
                # report P value
                outPvalue <- Rmimic::fuzzyP(pullvalue)
                if (outPvalue$modifier == "=") {
                  pullvalue <- outPvalue$report
                } else {
                  pullvalue <- sprintf('%s%s', outPvalue$modifier, outPvalue$report)
                }
                if (outPvalue$interpret <= studywiseAlpha) {
                  pullvalue <- sprintf('%s%s', pullvalue, "**")
                }
              }
              cat(sprintf("%-*s",sepgap[1,cC],pullvalue))
            } else {
              if ((!is.numeric(submatrix[1,sprintf("%s", vectnames[cC])])) & (!is.na(submatrix[1,sprintf("%s", vectnames[cC])]))) {
                cat(sprintf("%-*s",sepgap[1,cC],submatrix[1,sprintf("%s", vectnames[cC])]))
              } else {
                cat(sprintf("%-*s",sepgap[1,cC],"na"))
              }
            }
          }
          cat(sprintf("\n"))
        }
        if (cB < length(nongrouplabels)) {
          cat(sprintf("%s\n",paste(replicate(floor(spansize/3)-1, bigspancharacter), collapse = "")))
        }
      }
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      
      # output text writeups
      cat(sprintf("\n\n\n%s\n", "Interpretation with means, SD, and test statistics"))
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
      nongrouplabels <- unique(as.character(dataframeout$Variable))
      grouplabels <- unique(as.character(dataframeout$Comparison))
      ngroups <- length(grouplabels)
      for (cB in 1:length(nongrouplabels)) {
        submatrix <- dataframeout[which(dataframeout$Test == testlabels[cT]),] # subset for test
        submatrix <- submatrix[which(submatrix$Variable == nongrouplabels[cB]),] # subset for variable
        # Write Variable Label
        if (cB > 1) {
          cat(sprintf('\n'))
        }
        cat(sprintf("For %s:\n",nongrouplabels[cB]))
        for (cG in 1:nrow(submatrix)) {
          typewriter(submatrix$interpretation[cG], tabs=1, spaces=0, characters=spansize, indent="hanging")
          if (cG < nrow(submatrix)) {
            cat(sprintf('\n'))
          }
        } # end cG
      } # end cB
      cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    } # end cT
  } # end verbose
  
  # output to environmental variables
  res <- list()
  res$descriptives <- finalmasterdescriptives
  res$stats <- dataframeout
  return(res)
  
}
                             
                           