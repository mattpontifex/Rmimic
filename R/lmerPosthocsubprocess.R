#' lmerPosthocsubprocess
#'
#' @description Subfunction of Rmimic::lmerPosthoc. Performs posthoc tests.
#'
#' @param fit modelfit from the lmer function
#' @param dependentvariable text specifying the dependent variable label
#' @param subjectid text specifying the subject id label.
#' @param effectofinterest text specifying the effect label to decompose
#' @param within text list specifying any variables that should be considered as within subjects
#' @param between text list specifying the between subjects variables
#' @param covariates text list specifying the variables that should be considered as non interactive covariates. 
#' @param planned text list specifying any effect to show the post-hoc comparisons even if they are not significant.
#' @param df Parameter to indicate what degrees of freedom approximation should be used. Default is Kenward-Roger. Other options are Shattertwaite or Traditional.
#' @param confidenceinterval Decimal representation of confidence interval. Default 0.95.
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param posthoclimit integer specifying how many factors in an interaction should be broken down. Default is 6 indicating posthoc results will be provided for up to a 6 way interaction.
#' @param progressbar boolean parameter for if a progress bar should be shown indicating the status of the function. Note that multi way interactions may take a substantial period of time to fully decompose.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 22, 2025
#'
#' @importFrom progress progress_bar
#' @importFrom stringr str_split str_remove_all
#' @importFrom pkgcond suppress_conditions
#' @importFrom lmerTest lmer
#' @importFrom emmeans emmeans
#' @importFrom doBy summaryBy
#' @importFrom MBESS conf.limits.nct
#' @importFrom stats cor.test formula model.frame complete.cases
#' @importFrom WRS2 pbcor
#' @importFrom data.table data.table as.data.table
#' 
#' @export

lmerPosthocsubprocess <- function(results, effectofinterest, ...) {

  dots <- tryCatch({
    dots <- list(...) 
  }, error = function(e) {
    dots <- list() 
  })
  # pull info from results
  fit <- results$fit
  dependentvariable <- results$dependentvariable
  subjectid <- results$subjectid
  within <- results$within
  between <- results$between
  df <- results$df
  confidenceinterval <- results$confidenceinterval
  studywiseAlpha <- results$studywiseAlpha
  posthoclimit <- results$posthoclimit
  progressbar <- results$progressbar 
  covariates <- results$covariates 
  planned <- results$planned 
  bootstrap <- results$bootstrap
  
  debug <- FALSE
  if (debug) {
  
    dependentvariable <- 'Alertness'
    subjectid <- 'PartID'
    effectofinterest <- "Group:Time"
    within <- c('Time')
    df=NULL
    confidenceinterval=0.95
    studywiseAlpha=0.05
  
  }
  
  if (is.null(progressbar)) {
    progressbar <- FALSE
  }
  if (!is.null(df)) {
    if (toupper(df) == toupper("Kenward-Roger")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("KR")) {
      df = "Kenward-Roger"
    } else if (toupper(df) == toupper("Shattertwaite")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("S")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("satterthwaite")) {
      df = "Shattertwaite"
    } else if (toupper(df) == toupper("Traditional")) {
      df = "Traditional"
    }
  } else {
    df = "Kenward-Roger"
  }
  emmeansdf <- tolower(df)
  if (df == "Kenward-Roger") {
    emmeansdf <- 'kenward-roger'
  } else {
    emmeansdf <- 'satterthwaite'
  }
  
  # check that we are dealing with a lmer fit
  randomformula <- tryCatch({
    randomformula <- Reduce(paste, deparse(formula(fit, random.only = TRUE)))
  }, error = function(e) {
    randomformula <- NULL
  })
  if (!is.null(randomformula)) {
    # get data
    tempdbs <- tryCatch({
      tempdbs <- stats::model.frame(fit)
    }, error = function(e) {
      tempdbs <- NULL
    })
    
    factorsinvolved <- unlist(strsplit(as.character(effectofinterest),"[:]"))
    factorsinvolvedL <- length(factorsinvolved)
    
    # make sure this is something that we can do
    if (factorsinvolvedL < posthoclimit) {
        
      # see what we are working with
      factorlengthmatrix <- data.frame(matrix(NA,nrow=1, ncol=factorsinvolvedL))
      colnames(factorlengthmatrix) <- factorsinvolved
      for (cB in 1:factorsinvolvedL) {
        factorlengthmatrix[1,factorsinvolved[cB]] <- length(unique(unlist(as.character(tempdbs[,factorsinvolved[cB]]))))
      }
      
      res <- list()
      
      if (progressbar) {
        # establish progress
        formtext <- sprintf(" lmerPosthoc() decomposing %s [:bar] :percent eta: :eta", effectofinterest)
        #cat(sprintf('lmerPosthoc() beginning processing '))
        countticks <- 0
        stepcounts <- factorsinvolvedL
        pb <- progress::progress_bar$new(total = stepcounts,
                                         format = formtext,
                                         clear = TRUE, width= 120)
        pb$tick(0)
      }
      
      # loop through each factor
      for (currentFactorinvolved in 1:factorsinvolvedL) {
        currentfactor <- colnames(factorlengthmatrix)[currentFactorinvolved]
        currentfactorlevelsinvolved <- unique(unlist(as.character(tempdbs[,currentfactor])))
        otherfactorsinvolved <- factorsinvolved[which(factorsinvolved != currentfactor)]
        factortag <- currentfactor
        
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
          # ANOVA required
          for (cD in 1:length(currentfactorlevelsinvolved)) {
            # subset data
            subworkingdatabase <- tempdbs
            subworkingdatabase <- data.table::as.data.table(subworkingdatabase)
            subworkingdatabase <- subworkingdatabase[get(currentfactor[1]) == currentfactorlevelsinvolved[cD]]
            
            currentfactorlevelstring <- currentfactorlevelsinvolved[cD]
            # see if factor is just a number
            if (!grepl("\\D", currentfactorlevelstring)) {
              currentfactorlevelstring <- sprintf('%s%s', currentfactor[1], currentfactorlevelstring)
            }
            decompconst <- sprintf('When %s is %s', currentfactor[1], currentfactorlevelstring)
            
            # validate if random factors can be included
            rf_str <- Reduce(paste, deparse(stats::formula(fit, random.only = TRUE)))
            rf_terms <- stringr::str_split(stringr::str_remove_all(rf_str, " "), "~")[[1]][2]
            random_terms <- trimws(unlist(strsplit(rf_terms, "\\+")))
            keep_terms <- character(0)
            for (term in random_terms) {
              # Extract variable names
              vars <- unlist(regmatches(term, gregexpr("\\b\\w+\\b", term)))
              # Only keep variables that are columns in the data
              valid_vars <- intersect(vars, colnames(subworkingdatabase))
              # Vectorized unique count check
              uniq_counts <- sapply(valid_vars, function(v) data.table::uniqueN(subworkingdatabase[[v]]))
              # For all variables in term, must have at least 2 levels
              if (all(uniq_counts > 1)) {
                keep_terms <- c(keep_terms, term)
              } else {
                # Try replacing bad variables with 1 and clean up
                fixed_term <- term
                for (v in valid_vars[uniq_counts <= 1]) {
                  fixed_term <- stringr::str_replace_all(fixed_term, v, "1")
                  fixed_term <- stringr::str_replace_all(fixed_term, ":1|1:", "")
                }
                # If after replacement it's not trivial, keep it
                if (!fixed_term %in% keep_terms && fixed_term != "(1|1)" && fixed_term != "") {
                  keep_terms <- c(keep_terms, fixed_term)
                }
                # Otherwise, discard the term
              }
            }
            randomformula <- paste(keep_terms[keep_terms != ""], collapse = " + ")
            
            if (factorsinvolvedL > 1) {
              fixedformula <- paste(otherfactorsinvolved, collapse=sprintf("*")) # since this is an interaction decomposition
              if (length(otherfactorsinvolved) == 1) {
                #decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the effect of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), otherfactorsinvolved[1], currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), otherfactorsinvolved[1], currentfactorlevelstring)
              } else {
                #decomptext <- sprintf("Post-hoc decomposition of the %s interaction was conducted by examining the interaction of %s within each %s.", paste(factorsinvolved, collapse=sprintf(" \u00D7 ")), paste(otherfactorsinvolved, collapse=sprintf(" \u00D7 ")), currentfactor[1])
                factortag <- sprintf("%s_for_%s_within_%s", paste(factorsinvolved, collapse=sprintf("_by_")), 
                                     paste(otherfactorsinvolved, collapse=sprintf("_by_")), currentfactorlevelstring)
              }
            } else {
              fixedformula <- currentfactor[1]
              factortag <- currentfactor
            }
              
            # run the test on the factor of interest using the subsetted data
            smp <- as.data.frame(subworkingdatabase)
            subfit <- tryCatch({
              model_formula <- as.formula(sprintf('%s ~ %s + %s', dependentvariable[1], fixedformula, randomformula))
              subfit <- invisible(suppressWarnings(suppressMessages(pkgcond::suppress_conditions(lmerTest::lmer(formula = model_formula, data = smp)))))
            }, error = function(e) {
              cat(sprintf('\nThis error is usually the result of too many random factors specified for the interaction.\nUnable to run decomposition anova with the formula: %s ~ %s + %s\n', dependentvariable[1], fixedformula, randomformula))
              subfit <- NULL
            })
            if (!is.null(subfit)) {
              
              or_null <- function(x) if (length(x)) x else NULL
              normalize_vec <- function(x) if (is.null(x)) character(0) else as.character(x)
              keep_in_interest <- function(x, vars) {
                x    <- normalize_vec(x)
                vars <- normalize_vec(vars)
                x[x %in% vars]
              }
              smptempdbs <- stats::model.frame(subfit)
              varsofinterst <- colnames(smptempdbs)
              within <- or_null(keep_in_interest(within, varsofinterst))
              between <- or_null(keep_in_interest(between,    varsofinterst))
              covariates <- or_null(keep_in_interest(covariates, varsofinterst))
              
              # get table
              subresults <- invisible(suppressWarnings(suppressMessages(lmerEffects(subfit, dependentvariable=dependentvariable, subjectid = subjectid, within=within, df = df, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, suppresstext=FALSE, smp=smp))))
              # get posthoc
              subresults <- invisible(suppressWarnings(suppressMessages(lmerPosthoc(subresults, between=between, within=within, covariates=covariates, planned=planned, posthoccorrection='none', bootstrap=bootstrap, confidenceinterval=confidenceinterval, studywiseAlpha=studywiseAlpha, posthoclimit=posthoclimit-1, calltype='subprocess', progressbar=progressbar))))
              
              # need to output subresults to res
              if (length(names(subresults)) > 0) {
                res[[sprintf("ANOVA_%s", factortag)]] <- subresults
              }
            } # submodel
          } # each factor level involved
        } else {
          # just normal results
          outtable <- lmerPosthocsubprocessContrasts(fit, emmeansdf, effectofinterest, currentfactor, factorsinvolved, otherfactorsinvolved, dependentvariable, decomptext, subjectid, confidenceinterval, studywiseAlpha, within, factortag)

          # check to see if bootstrapping should be performed
          if (!is.null(bootstrap)) {
            if (!isFALSE(bootstrap)) {
              contrastin <- list()
              contrastin$emmeansdf <- emmeansdf
              contrastin$effectofinterest <- effectofinterest
              contrastin$currentfactor <- currentfactor
              contrastin$factorsinvolved <- factorsinvolved
              contrastin$otherfactorsinvolved <- otherfactorsinvolved
              contrastin$decomptext <- decomptext
              contrastin$factortag <- factortag
              contrastin$outtable <- outtable
              outtable <- lmerEffectsBootstrapContrast(results, contrastin, bootstrap)
            }
          }
          # need to output table to res
          if (nrow(outtable) > 0) {
            res[[sprintf("Posthoc_%s", factortag)]] <- outtable
          }
        } # no anova required
        
        if (progressbar) {
          # update progress
          #  cat(sprintf('.'))
          pb$tick()
          countticks <- countticks + 1
        }
      }
      if (progressbar) {
        if (countticks < stepcounts) {
          pb$tick(stepcounts)
        }
      }
    } # posthoc limit
  }
  
  return(res)
}