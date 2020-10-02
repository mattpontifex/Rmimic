#' pop_RmimicMediation
#'
#' @description Popup window to compute mediation analysis with confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, October 1, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer HTML reactive
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' @importFrom stats glm lm binomial
#' @importFrom mediation mediate
#' 
#'
#' @export

pop_RmimicMediation <- function() {
  
  # Get the document context.
  #context <- rstudioapi::getActiveDocumentContext()
  
  # Extract current data from environment
  environelements <- ls(envir=.GlobalEnv) # global elements
  info_dfs <- ""
  if (length(environelements) > 0) {
    environelementsidx <- which(sapply(environelements, function(x) is.data.frame(get(x)))) # only dataframes
    environelements <- environelements[environelementsidx] 
    
    info_dfs <- lapply(
      X = environelements,
      FUN = function(x) {
        tmp <- get(x, envir = .GlobalEnv)
        sprintf("%d obs. of  %d variables", nrow(tmp), ncol(tmp))
      }
    )
    info_dfs <- unlist(info_dfs)
  }
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Rmimic: Compute Mediation", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
    miniUI::miniContentPanel(
      padding = c(30, 20, 20, 20),
      shinyWidgets::pickerInput(
        inputId = "select_dataframe",
        label = "Select your working data.frame:",
        choices = environelements, width = "100%",
        options = list(title = "List of data.frame..."),
        choicesOpt = list(subtext = info_dfs),
        selected = NULL
      ),
      
      shiny::uiOutput("ui1") # This outputs the dynamic UI component
    ),
    
    shinyWidgets::actionBttn(
      inputId = "generatecode",
      label = "Generate Code Only",
      style = "simple", 
      color = "royal",
      block=TRUE,
      size="sm"
    ),
    
    shinyWidgets::actionBttn(
      inputId = "done",
      label = "Compute (this may take a moment)",
      style = "simple", 
      color = "primary",
      block=TRUE,
      size="md"
    )
    
    
  )
  
  server <- function(input, output) {
    
    output$ui1 <- shiny::renderUI({
      if (input$select_dataframe == "") {
        return()
      } else {
        
        # Extract current data from environment
        workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
        workingdatavariables <- names(workingdata)
        names(workingdatavariables) <- workingdatavariables
        workingdatavariablestypes <- c()
        
        for (cC in 1:length(workingdatavariables)) {
          testvect <- unlist(workingdata[,cC])
          testvectout <- NA
          if (!is.null(testvect)) {
            if (inherits(testvect, c("Date", "POSIXct", "POSIXlt"))) {
              testvectout <- "time/date"
            } else if (inherits(testvect, c("logical", "factor", "AsIs"))) {
              testvectout <- "logical"
            } else if (inherits(testvect, c("character"))) {
              testvectout <- "character"
            } else if (inherits(testvect, c("numeric", "integer", "double"))) {
              testvectout <- "numeric"
            }
          }
          workingdatavariablestypes <- c(workingdatavariablestypes, sprintf('    type: %s', testvectout))
        }
        
        workingdatavariables <- workingdatavariables
        workingdatavariablestypes <- workingdatavariablestypes
        
        shiny::wellPanel(
          shiny::tagList(
            
            shinyWidgets::pickerInput(
              inputId = "select_DV",
              label = "Select the Outcome variable (dependent variable):",
              choices = workingdatavariables, width = "100%",
              multiple = FALSE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::pickerInput(
              inputId = "select_IV",
              label = 'Select the Predictor variable (independent variable)',
              choices = workingdatavariables, width = "100%",
              multiple = FALSE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::pickerInput(
              inputId = "select_MV",
              label = 'Select the Mediating variable',
              choices = workingdatavariables, width = "100%",
              multiple = FALSE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::pickerInput(
              inputId = "select_CV",
              label = 'Select any other variables to account for',
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_modelstyle",
              label = "What type of regression should be run?",
              choices = c("Linear Regression","Logistic Regression"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "100%",
            ),
            
            shinyWidgets::sliderTextInput(
              inputId = "select_simnum",
              label = "How many Monte Carlo draws should be used?", 
              choices = seq(from = 500, to = 5000, by = 25),
              selected = 1000,
              width = "80%"
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_bootstyle",
              label = "What type of bootstrapping should be run?",
              choices = c("Quasi-Bayesian","Nonparametric"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "100%",
            ),
            
            
          )
          
        )
      }
    })
    
    toListen <- shiny::reactive({
      list(input$done,input$generatecode)
    })
    
    shiny::observeEvent(toListen(), {
      if ((input$generatecode != 0) | (input$done != 0)) {
        
        # Check that selections are made
        if (!is.null(input$select_dataframe)) {
          if (!is.null(input$select_DV)) {
            if (!is.null(input$select_IV)) {
              if (!is.null(input$select_MV)) {
            
                listofcalls <- c()
                tmppref <- ''
                tempsuff <- ''
                if (input$select_modelstyle == 'Linear Regression') {
                  tmppref <- 'fitM <- stats::lm('
                  tmppref2 <- 'fitY <- stats::lm('
                  tempsuff <- ''
                } else {
                  tmppref <- 'fitM <- stats::glm('
                  tmppref2 <- 'fitY <- stats::glm('
                  tempsuff <- 'family=stats::binomial(link = "logit"), '
                }
                tempsuff <- sprintf('%sdata=%s)', tempsuff, input$select_dataframe)
              
                listofvariables <- input$select_IV[1]
                if (length(input$select_CV) > 0) {
                  listofvariables <- c(input$select_IV[1], input$select_CV)
                }
                listofvariables <- paste(listofvariables, collapse=" + ")
                tmpform <- sprintf('%s ~ %s', input$select_MV[1], listofvariables) 
                tmpcall <- sprintf('%s%s, data=%s)',tmppref,tmpform,input$select_dataframe) 
                listofcalls <- c(listofcalls, tmpcall)
                
                listofvariables <- c(input$select_IV[1], input$select_MV[1])
                if (length(input$select_CV) > 0) {
                  listofvariables <- c(input$select_IV[1], input$select_MV[1], input$select_CV)
                }
                listofvariables <- paste(listofvariables, collapse=" + ")
                tmpform <- sprintf('%s ~ %s', input$select_DV[1], listofvariables) 
                tmpcall <- sprintf('%s%s, %s',tmppref2,tmpform,tempsuff) 
                listofcalls <- c(listofcalls, tmpcall)
                
                tmpcall <- sprintf('fitMed <- mediation::mediate(fitM, fitY')
                tmpcall <- sprintf('%s, treat=%s', tmpcall, sprintf("'%s'",input$select_IV[1]))
                tmpcall <- sprintf('%s, mediator=%s', tmpcall, sprintf("'%s'",input$select_MV[1]))
                if (input$select_bootstyle == "Nonparametric") {
                  tmpcall <- sprintf('%s, boot=TRUE', tmpcall)
                } else {
                  tmpcall <- sprintf('%s, boot=FALSE', tmpcall)
                }
                tmpcall <- sprintf('%s, sims=%d', tmpcall,input$select_simnum)
                tmpcall <- sprintf('%s, conf.level=0.95)', tmpcall)
                listofcalls <- c(listofcalls, tmpcall)
                
                tmpcall <- sprintf('res <- Rmimic::mediate2text(fitMed, studywiseAlpha=0.05)')
                listofcalls <- c(listofcalls, tmpcall)
                
                # execute call
                codelevel <- 0 
                if (input$done) {
                  boolattempt <- FALSE
                  for (cR in 1:length(listofcalls)) {
                    tmpcall <- listofcalls[cR]
                    boolattempt <- tryCatch({
                      eval(parse(text=tmpcall))
                      boolattempt <- TRUE}
                    )
                  }
                  if (boolattempt == FALSE) {
                    Rmimic::typewriter('Uh oh.. Something went wrong. But the syntax for the function is provided below.', tabs=0, spaces=0, characters=80, indent='hanging')
                  }
                  
                  # output calls
                  Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=80, indent='hanging')
                  codelevel <- 1
                }
                for (cR in 1:length(listofcalls)) {
                  tmpcall <- listofcalls[cR]
                  Rmimic::typewriter(tmpcall, tabs=codelevel, spaces=0, characters=80, indent='hanging')
                }
              }
            }
          } # DV select
        } # dataframe select
        
        invisible(shiny::stopApp())
      }
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}
