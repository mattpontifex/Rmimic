#' mediate2text
#'
#' @description Output mediate results in a more intelligible format.
#'
#' @param fit mediate model fit
#' @param studywiseAlpha Decimal representation of alpha level. Default 0.05.
#' @param verbose Boolean operator for if interpretations of the statistics should be printed. Default is TRUE.
#' 
#' @return A list with:
#' \item{fit}{mediate object}
#' \item{summary}{summary object}
#' \item{paths}{table for path effects}
#' \item{interpretation}{text interpretation}
#' \item{propmediated}{proportion of the effect mediated}
#' 
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 1, 2020
#' 
#' @examples
#' 
#' 
#' @importFrom stats confint 
#'
#' @export

mediate2text <- function(fit, studywiseAlpha=0.05, verbose=TRUE) {
  
  spancharacter <- "-"
  spansize <- 95
  operatingsystem <- Sys.info()['sysname']
  if (operatingsystem == "Windows") {
    spancharacter <- "_"
  } else if (operatingsystem == "Darwin") {
    spancharacter <- "-"
  }  
  
  confidenceinterval <- fit$conf.level
  Mediatorlabel <- fit$mediator
  IVlabel <- fit$treat
  fitMS <- summary(fit)
  DVlabel <- toString(fitMS[["model.y"]][["terms"]][[2]])
  
  result <- list()
  result$fit <- fit
  result$summary <- fitMS
  result$paths <- data.frame(matrix(NA, nrow=4, ncol = 5))
  names(result$paths) <- c('path', 'estimate', 'conf.int.lower', 'conf.int.upper', 'p.value')
  
  outstring <- sprintf('Direct Path:')
  outstring <- sprintf('%s %s ---> %s', outstring, IVlabel, DVlabel)
  result$paths[1,1] <- outstring
  #(taken from fitMS Total Effect)
  result$paths[1,2] <- fitMS[["tau.coef"]]
  result$paths[1,3] <- fitMS[["tau.ci"]][[1]]
  result$paths[1,4] <- fitMS[["tau.ci"]][[2]]
  result$paths[1,5] <- fitMS[["tau.p"]]
  
  outstring <- sprintf('Mediated Path A: ')
  outstring <- sprintf('%s %s ---> %s', outstring, IVlabel, Mediatorlabel)
  result$paths[2,1] <- outstring
  result$paths[2,2] <- fit$model.m[["coefficients"]][[IVlabel]]
  
  tempA <- stats::confint(fit$model.m, level=confidenceinterval)
  tempAi <- which(rownames(tempA) == IVlabel)
  result$paths[2,3] <- tempA[tempAi,1]
  result$paths[2,4] <- tempA[tempAi,2]
  msA <- summary(fit$model.m)
  msAf <- data.frame(msA$coefficients)
  tempAi <- which(rownames(msAf) == IVlabel)
  result$paths[2,5] <- msAf[tempAi,4]
  
  outstring <- sprintf('Mediated Path B: ')
  outstring <- sprintf('%s %s ---> %s', outstring, Mediatorlabel, DVlabel)
  result$paths[3,1] <- outstring
  #(taken from fitMS ACME)
  result$paths[3,2] <- fitMS[["d0"]]
  result$paths[3,3] <- fitMS[["d0.ci"]][[1]]
  result$paths[3,4] <- fitMS[["d0.ci"]][[2]]
  result$paths[3,5] <- fitMS[["d0.p"]]
  
  outstring <- sprintf('Mediated Path C: ')
  outstring <- sprintf('%s %s ---> %s', outstring, IVlabel, DVlabel)
  result$paths[4,1] <- outstring
  #(taken from fitMS ADE)
  result$paths[4,2] <- fitMS[["z0"]]
  result$paths[4,3] <- fitMS[["z0.ci"]][[1]]
  result$paths[4,4] <- fitMS[["z0.ci"]][[2]]
  result$paths[4,5] <- fitMS[["z0.p"]]
  
  templist <- stringr::str_split(toString(fit$model.y$terms), ",")[[1]]
  templist <- stringr::str_trim(templist)[3]
  templist <- stringr::str_split(templist, " ")[[1]]
  templist <- templist[which(templist != "+")]
  templist <- templist[which(templist != IVlabel)]
  templist <- templist[which(templist != Mediatorlabel)]
  covariates <- ""
  covariateL <- 0
  if (length(templist) > 0) {
    covariates <- templist
    covariateL <- length(templist)
  }
  
  if (verbose == TRUE) {
    
    temptext <- "Mediation Analysis"
    temptextspan <- floor(nchar(temptext)/2)
    pagespan <- floor(spansize/2)
    cat(sprintf("\n"))
    Rmimic::typewriter(temptext, tabs=0, spaces=(pagespan-temptextspan), characters=floor(spansize*.9))
    
    outstring <- "Mediation analyses were conducted using the"
    outstring <- sprintf('%s mediation (Tingley,Yamamoto, Hirose, Keele, & Imai, 2014)', outstring)
    outstring <- sprintf('%s and Rmimic (Pontifex, %s) packages', outstring, strsplit(as.character(utils::packageDate("Rmimic")),"-")[[1]][1])
    rvers <- unlist(strsplit(R.version.string, " "))
    rvers <- paste(rvers[1:length(rvers)-1], collapse=" ")
    outstring <- sprintf('%s in %s.', outstring, rvers)
    
    outstring <- sprintf('%s Unstandardized indirect effects were computed using', outstring)
    if (fit$sims > 0) {
      outstring <- sprintf('%s %d', outstring, fit$sims)
    }
    if (fit$boot == TRUE) {
      outstring <- sprintf('%s nonparametric bootstraped', outstring)
    } else {
      outstring <- sprintf('%s quasi-Bayesian approximation based', outstring)
    }
    outstring <- sprintf('%s samples', outstring)
    if (fit$boot == TRUE) {
      if (fit$robustSE == TRUE) {
        outstring <- sprintf('%s with heteroskedasticity-consistent standard errors', outstring)
      }
    }
    outstring <- sprintf('%s.', outstring)
    
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
    cat(sprintf("\n"))
    
    outstring <- sprintf('%s', result$paths[1,1])
    Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
    
    outstring <- sprintf('Estimate = %0.2f', round(result$paths[1,2],2))
    outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(result$paths[1,3],2), round(result$paths[1,4],2))
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$paths[1,5]))
    outstring <- sprintf('%s, p %s %s', outstring, outPvalue$modifier, outPvalue$report)
    if (outPvalue$interpret <= studywiseAlpha) {
      outstring <- sprintf('%s%s', outstring, "**")
    }
    Rmimic::typewriter(outstring, tabs=4, spaces=0, characters=floor(spansize*.9))
    cat(sprintf("\n"))
    
    outstring <- sprintf('Mediated Path:')
    Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
    outstring <- sprintf('Path A: ')
    outstring <- sprintf('%s %s ---> %s', outstring, IVlabel, Mediatorlabel)
    Rmimic::typewriter(outstring, tabs=2, spaces=0, characters=floor(spansize*.9))

    outstring <- sprintf('Estimate = %0.2f', round(result$paths[2,2],2))
    outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(result$paths[2,3],2), round(result$paths[2,4],2))
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$paths[2,5]))
    outstring <- sprintf('%s, p %s %s', outstring, outPvalue$modifier, outPvalue$report)
    if (outPvalue$interpret <= studywiseAlpha) {
      outstring <- sprintf('%s%s', outstring, "**")
    }
    Rmimic::typewriter(outstring, tabs=4, spaces=0, characters=floor(spansize*.9))
    cat(sprintf("\n"))
    
    outstring <- sprintf('Path B: ')
    outstring <- sprintf('%s %s ---> %s', outstring, Mediatorlabel, DVlabel)
    Rmimic::typewriter(outstring, tabs=2, spaces=0, characters=floor(spansize*.9))
    outstring <- sprintf('Estimate = %0.2f', round(result$paths[3,2],2))
    outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(result$paths[3,3],2), round(result$paths[3,4],2))
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$paths[3,5]))
    outstring <- sprintf('%s, p %s %s', outstring, outPvalue$modifier, outPvalue$report)
    if (outPvalue$interpret <= studywiseAlpha) {
      outstring <- sprintf('%s%s', outstring, "**")
    }
    Rmimic::typewriter(outstring, tabs=4, spaces=0, characters=floor(spansize*.9))
    cat(sprintf("\n"))
    
    outstring <- sprintf('Path C: ')
    outstring <- sprintf('%s %s ---> %s', outstring, IVlabel, DVlabel)
    Rmimic::typewriter(outstring, tabs=2, spaces=0, characters=floor(spansize*.9))
    outstring <- sprintf('Estimate = %0.2f', round(result$paths[4,2],2))
    outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(result$paths[4,3],2), round(result$paths[4,4],2))
    outPvalue <- Rmimic::fuzzyP(as.numeric(result$paths[4,5]))
    outstring <- sprintf('%s, p %s %s', outstring, outPvalue$modifier, outPvalue$report)
    if (outPvalue$interpret <= studywiseAlpha) {
      outstring <- sprintf('%s%s', outstring, "**")
    }
    Rmimic::typewriter(outstring, tabs=4, spaces=0, characters=floor(spansize*.9))
    
    if (covariateL > 0) {
      cat(sprintf("\n"))
      outstring <- sprintf('Accounting for the influence of: %s', paste(covariates, collapse=", "))
      Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
    }
    
    outstring <- "Interpretation with Test Statistics"
    cat(sprintf("\n"))
    Rmimic::typewriter(outstring, tabs=0, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
    
  }
      
  outstring <- sprintf('The relationship between %s and %s was', IVlabel, DVlabel) 
  outPvalueACME <- Rmimic::fuzzyP(as.numeric(fitMS[["d0.p"]]))
  outPvalueADE <- Rmimic::fuzzyP(as.numeric(fitMS[["z0.p"]]))
  if (outPvalueACME$interpret <= studywiseAlpha) {
    if (outPvalueADE$interpret <= studywiseAlpha) {
      outstring <- sprintf('%s partially', outstring)
    } else {
      outstring <- sprintf('%s fully', outstring)
    }
  } else {
    outstring <- sprintf('%s not', outstring) 
  }
  outstring <- sprintf('%s mediated by %s', outstring, Mediatorlabel) 
  
  outstring <- sprintf('%s (Proportion Mediated = %0.1f%%', outstring, round(fitMS[["n0"]],3)*100)
  
  outstring <- sprintf('%s; Average Causal Mediation Effect = %0.2f', outstring, round(fitMS[["d0"]], 2))
  outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(fitMS[["d0.ci"]][[1]],2), round(fitMS[["d0.ci"]][[2]],2))
  outstring <- sprintf('%s, p %s %s', outstring, outPvalueACME$modifier, outPvalueACME$report)
  
  outstring <- sprintf('%s; Average Direct Effect = %0.2f', outstring, round(fitMS[["z0"]], 2))
  outstring <- sprintf('%s [%2.0f%% CI: %0.2f to %0.2f]', outstring, floor(confidenceinterval*100), round(fitMS[["z0.ci"]][[1]],2), round(fitMS[["z0.ci"]][[2]],2))
  outstring <- sprintf('%s, p %s %s)', outstring, outPvalueADE$modifier, outPvalueADE$report)
  
  if (covariateL > 0) {
    cat(sprintf("\n"))
    outstring <- sprintf('%s, after controlling for the effects of', outstring)
    
    if (covariateL == 1) {
      outstring <- sprintf('%s %s', outstring, covariates[1])
    } else if (covariateL == 2) {
      outstring <- sprintf('%s %s and %s', outstring, covariates[1], covariates[2])
    } else {
      for (cInc in 1:covariateL) {
        outstring <- sprintf('%s %s', outstring, covariates[cInc])
        if (cInc < covariateL) {
          outstring <- sprintf('%s,', outstring)
        }
        if (cInc == (covariateL-1)) {
          outstring <- sprintf('%s and', outstring)
        }
      }
    }
    
  }
  outstring <- sprintf('%s.', outstring)
  
  result$interpretation <- outstring
  result$propmediated <- fitMS[["n0"]]
  
  if (verbose == TRUE) {
    Rmimic::typewriter(outstring, tabs=1, spaces=0, characters=floor(spansize*.9))
    rm(outstring)
    cat(sprintf("\n"))
    cat(sprintf("%s\n",paste(replicate(spansize, spancharacter), collapse = "")))
  }
  
  return(result)
}
