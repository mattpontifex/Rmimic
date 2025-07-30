#' lmerPosthocCorrection
#'
#' @description Correction of posthoc tests. Subfunction of lmerPosthoc. Corrections are made in place and returned. 
#'
#' @param results list containing data frames for stats, randomstats, and rsquared values
#' @param method text indicating the correction to apply. Options are Bonferroni, Sidak, Holm-Bonferroni, or False Discovery Rate Control (default)
#' @param studywiseAlpha decimal indicating the alpha level
#' @param FDRC decimal specifying the false discovery rate control. Default is 0.05
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 1, 2025
#' 
#' @examples
#'
#'     fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = Rmimic::alertness)
#'     results <- Rmimic::lmerEffects(fit, dependentvariable = "Alertness", subjectid = "PartID", df = "Kenward-Roger")
#'     results <- Rmimic::lmerPosthoc(results, between=c('Group'), within=c('Time'),
#'                covariates=NULL, planned=c('Group'), posthoccorrection="none", progressbar=TRUE)
#'     results <- Rmimic::lmerPosthocCorrection(results, method="False Discovery Rate Control", studywiseAlpha=0.05, FDRC=0.05)      
#'                
#' @export

lmerPosthocCorrection <- function(results, method=NULL, studywiseAlpha=NULL, FDRC=0.05) {
  
  if (is.null(studywiseAlpha)) {
    studywiseAlpha <- tryCatch({
      studywiseAlpha <- results$studywiseAlpha
    }, error = function(e) {
      studywiseAlpha <- NULL
    })
    if (is.null(studywiseAlpha)) {
      studywiseAlpha <- 0.05
      results$studywiseAlpha <- studywiseAlpha # store it
    }
  }
  
  if (!is.null(method)) {
    if (toupper(method) == toupper("Bonferroni")) {
      method = "Bonferroni"
    } else if (toupper(method) == toupper("Sidak")) {
      method = "Sidak"
    } else if (toupper(method) == toupper("Holm-Bonferroni")) {
      method = "Holm-Bonferroni"
    } else if (toupper(method) == toupper("Holm")) {
      method = "Holm-Bonferroni"
    } else if (toupper(method) == toupper("None")) {
      method = "None"
    } else if (toupper(method) == toupper("False Discovery Rate Control")) {
      method = "False Discovery Rate Control"
    }
  } else {
    method = "False Discovery Rate Control"
  }
  
  #"Bonferroni", "Sidak","Holm-Bonferroni" are addressed within the t-test
  if ((toupper(method) == toupper("Bonferroni")) | (toupper(method) == toupper("Sidak")) | (toupper(method) == toupper("Holm-Bonferroni"))) {
    
    if ('posthoc' %in% names(results)) {
      uniquenames <- names(results$posthoc)
      for (cUN in 1:length(uniquenames)) {
  
        #textcall <- sprintf("tempelement <- results$posthoc$%s", uniquenames[cUN])
        #eval(parse(text=textcall))
        tempelement <- results$posthoc[[uniquenames[cUN]]]
        
        if (is.data.frame(tempelement)) {
          # element is a table
          
          if ((toupper(method) == toupper("Bonferroni")) | (toupper(method) == toupper("Sidak"))) {
            critp <- studywiseAlpha
            criticalphrase <- ""
            if (toupper(method) == toupper("Bonferroni")) {
              critp <- (studywiseAlpha/nrow(tempelement))
              criticalphrase <- sprintf("(Bonferroni critical alpha = %.4f)", round(critp,4))
            } else if (toupper(method) == toupper("Sidak")) {
              critp <- (1-(1-studywiseAlpha)^(1/nrow(tempelement)))
              criticalphrase <- sprintf("(Sidak critical alpha = %.4f)", round(critp,4))
            }
            critp <- round(round(round(critp, digits=7), digits=6), digits=5)
            for (cR in 1:nrow(tempelement)) {
              if (Rmimic::fuzzyP(tempelement$p.value[cR], studywiseAlpha)$significance) {
                if (!(Rmimic::fuzzyP(tempelement$p.value[cR], critp)$significance)) {
                  tempelement$textoutput[cR] <- sprintf('%s <span class="posthoccorrection">However, that difference did not remain significant following correction for multiple comparisons, %s.</span>', tempelement$textoutput[cR], criticalphrase)
                }
              }
            }
          }
          
          if (toupper(method) == toupper("Holm-Bonferroni")) {
            #Holm, S. 1979. A simple sequential rejective multiple test procedure. Scandinavian Journal of Statistics 6:65-70
            temp <- data.frame(matrix(NA, nrow=nrow(tempelement), ncol=2))
            colnames(temp) <- c("p.value", "location")
            for (cR in 1:nrow(tempelement)) {
              outPvalue <- Rmimic::fuzzyP(tempelement$p.value[cR], studywiseAlpha)
              if (outPvalue$significance) {
                temp[cR,1] <- outPvalue$exact
                temp[cR,2] <- cR
              }
            }
            temp <- temp[order(-temp$p.value),]
            ncomp <- length(which(temp$p.value <= as.numeric(studywiseAlpha)))
            
            # Loop through P values
            rank <- 1
            while (rank <= nrow(temp)) {
              if (!is.na(temp$p.value[rank])) {
                if (Rmimic::fuzzyP(temp$p.value[rank], studywiseAlpha)$significance) {
                  temppval <- (studywiseAlpha/(ncomp - rank + 1))
                  temppval <- round(round(round(temppval, digits=7), digits=6), digits=5)
                  if (!(Rmimic::fuzzyP(temp$p.value[rank], temppval)$significance)) {
                    # P value is no longer considered significant
                    criticalphrase <- sprintf("(Holm-Bonferroni critical alpha = %.4f)", round(temppval, digits=4))
                    tempelement$textoutput[temp$location[rank]] <- sprintf('%s <span class="posthoccorrection">However, that difference did not remain significant following correction for multiple comparisons, %s.</span>', tempelement$textoutput[temp$location[rank]], criticalphrase)
                  }
                } else {
                  # P value is not significant
                  rank <- rank + 1
                }
              }
              rank <- rank + 1
            }
          }
          
        } else {
          # element is a structure
          # recursive call
          tempelement <- lmerPosthocCorrection(tempelement, method=method, studywiseAlpha=studywiseAlpha)
        }
        
        # put it back
        #textcall <- sprintf("results$posthoc$%s <- tempelement", uniquenames[cUN])
        #eval(parse(text=textcall))
        results$posthoc[[uniquenames[cUN]]] <- tempelement
        
      }
    }
  }
  
  if (toupper(method) == toupper("False Discovery Rate Control")) {
    # Glickman, M. E., Rao, S. R., Schultz, M. R. (2014). False discovery rate control is a recommended alternative to Bonferroni-type adjustments in health studies. Journal of Clinical Epidemiology, 67, 850-857.
    
    if ((!is.numeric(FDRC)) | ((FDRC < 0) | (FDRC > 0.99))) {
      FDRC <- 0.05
    }
    # get posthoclinkmatrix
    results <- lmerPosthocCorrectionsubprocess(results)
    
    if (nrow(results$posthoclinkmatrix) > 0) {
      temp <- results$posthoclinkmatrix[order(results$posthoclinkmatrix$p.value),]
      ncomp <- nrow(temp)
      # Loop through P values
      for (rank in 1:nrow(temp)) {
        if (temp$significant[rank]) {
          temppval <- (FDRC*(rank/ncomp))
          temppval <- round(round(round(temppval, digits=7), digits=6), digits=5)
          if (!(Rmimic::fuzzyP(as.numeric(temp$p.value[rank]), temppval)$significance)) {
            # P value is no longer considered significant
            criticalphrase <- sprintf("(Benjamini-Hochberg critical alpha = %.4f)", round(temppval, digits=4))
            # get text
            templocation <- temp$location[rank]
            temprow <- temp$row[rank]
            #funcal <- sprintf('tempinterlocation <- results$posthoc$%s$textoutput[%d]',templocation,temprow)
            #suppressWarnings(eval(parse(text=funcal)))
            tempinterlocation <- results$posthoc[[sprintf('%s$textoutput[%d]',templocation,temprow)]]
            
            # ammend text with warning
            tempinterlocation <- sprintf('%s <span class="posthoccorrection">However, that difference did not remain significant following false discovery rate control %s.</span>', tempinterlocation, criticalphrase)
            # put it back
            #funcal <- sprintf('results$posthoc$%s$textoutput[%d] <- tempinterlocation',templocation,temprow)
            #suppressWarnings(eval(parse(text=funcal)))
            results$posthoc[[sprintf('%s$textoutput[%d]',templocation,temprow)]] <- tempinterlocation
          }
        }
      }
    }
  }
  
  return(results)
}