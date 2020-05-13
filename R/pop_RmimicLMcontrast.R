#' pop_RmimicLMcontrast
#'
#' @description Popup window to compute SPSS style regression analysis with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 5, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#'
#' @export

pop_RmimicLMcontrast <- function() {
  
  # Get the document context.
  #context <- rstudioapi::getActiveDocumentContext()
  
  # Extract current data from environment
  environelements <- ls(envir=.GlobalEnv) # global elements
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
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Rmimic: Compute Regression", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
      inputId = "done",
      label = "Compute",
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
        workingdatavariables <- c('< constant >', workingdatavariables)
        workingdatavariablestypes <- c('    type: constant', workingdatavariablestypes)
        
        model1contrasttext <- HTML("Model 1: Select the explanatory variables (independent variables) for the base model:<br><small>
                                   This can be a constant or control variables.</small>")
        
        shiny::wellPanel(
          shiny::tagList(
            
            shinyWidgets::pickerInput(
              inputId = "select_DV",
              label = "Select the variable you are trying to predict (dependent variable):",
              choices = workingdatavariables, width = "100%",
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
            
            shinyWidgets::pickerInput(
              inputId = "select_M1variables",
              label = model1contrasttext,
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_M1style",
              label = "How should Model 1 be run?",
              choices = c("All Variables Entered","Forward Stepwise","Backward Stepwise"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "100%",
            ),
            
            
            shinyWidgets::pickerInput(
              inputId = "select_M2variables",
              label = "Model 2: Select the explanatory variables (independent variables) for the model of interest:",
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_M2style",
              label = "How should Model 2 be run?",
              choices = c("All Variables Entered","Forward Stepwise","Backward Stepwise"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "100%",
            ),
            
          )
          
        )
      }
    })
    
    shiny::observeEvent(input$done, {
      
      # Check that selections are made
      if (!is.null(input$select_dataframe)) {
        if (!is.null(input$select_DV)) {
          # user has chosen a DV
          if ((!is.null(input$select_M1variables)) | (!is.null(input$select_M2variables))) {
            
            print("line171")
            M1variables <- input$select_M1variables
            if (length(which(M1variables == '< constant >')) > 0) {
              M1variables[which(M1variables == '< constant >')] <- '1'
            }
            M2variables <- input$select_M2variables
            if (length(which(M2variables == '< constant >')) > 0) {
              M2variables[which(M2variables == '< constant >')] <- '1'
            }
            if (length(M1variables) == 0) {
              # if no model 1 was entered, put constant in
              M1variables <- '1'
            }
            if (length(M2variables) == 0) {
              # if no model 2 was entered, put constant in model 1 and shift
              M2variables <- M1variables
              M1variables <- '1'
            }
            
            listofcalls <- c()
            tmppref <- ''
            tempsuff <- ''
            if (input$select_modelstyle == 'Linear Regression') {
              tmppref <- 'basefit <- stats::lm('
              tmppref2 <- 'fit <- stats::lm('
              tempsuff <- ''
            } else {
              tmppref <- 'basefit <- stats::glm('
              tmppref2 <- 'fit <- stats::glm('
              tempsuff <- 'family=binomial(link = "logit"), '
            }
            tempsuff <- sprintf('%sdata=%s)', tempsuff, input$select_dataframe)
            
            # Model 1
            tmpform <- sprintf('%s ~ %s', input$select_DV[1], paste(M1variables, collapse=" + "))
            
            if (input$select_M1style == "All Variables Entered") {
              tmpcall <- sprintf('%s%s, %s',tmppref,tmpform,tempsuff)
              listofcalls <- c(listofcalls, tmpcall)
            } else if (input$select_M1style == "Forward Stepwise") {
              tmpformsimple <- sprintf('%s ~ 1', input$select_DV)
              tmpcallsimple <- sprintf('%s%s, %s',tmppref,tmpformsimple,tempsuff)
              listofcalls <- c(listofcalls, tmpcallsimple)
              tmpcallalt <- sprintf('basefit <- stats::update(basefit,formula=MASS::stepAIC(basefit, direction="forward", scope=%s, trace=TRUE)$terms)', tmpform)
              listofcalls <- c(listofcalls, tmpcallalt)
            } else if (input$select_M1style == "Backward Stepwise") {
              tmpcall <- sprintf('%s%s, %s',tmppref,tmpform,tempsuff)
              listofcalls <- c(listofcalls, tmpcall)
              tmpcallalt <- sprintf('basefit <- stats::update(basefit,formula=MASS::stepAIC(basefit, direction="backward", trace=TRUE)$terms)')
              listofcalls <- c(listofcalls, tmpcallalt)
            }
            
            # Model 2
            tmpform <- sprintf('%s ~ %s', input$select_DV[1], paste(M2variables, collapse=" + "))
            if (input$select_M2style == "All Variables Entered") {
              tmpcall <- sprintf('%s%s, %s',tmppref2,tmpform,tempsuff)
              listofcalls <- c(listofcalls, tmpcall)
            } else if (input$select_M2style == "Forward Stepwise") {
              tmpformsimple <- sprintf('%s ~ 1', input$select_DV)
              tmpcallsimple <- sprintf('%s%s, %s',tmppref2,tmpformsimple,tempsuff)
              listofcalls <- c(listofcalls, tmpcallsimple)
              tmpcallalt <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, direction="forward", scope=%s, trace=TRUE)$terms)', tmpform)
              listofcalls <- c(listofcalls, tmpcallalt)
            } else if (input$select_M2style == "Backward Stepwise") {
              tmpcall <- sprintf('%s%s, %s',tmppref2,tmpform,tempsuff)
              listofcalls <- c(listofcalls, tmpcall)
              tmpcallalt <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, direction="backward", trace=TRUE)$terms)')
              listofcalls <- c(listofcalls, tmpcallalt)
            }
            
            tmpcall <- 'regresult <- Rmimic::RmimicLMcontrast(basefit, fit, studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE)'
            listofcalls <- c(listofcalls, tmpcall)
            
            # execute call
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
            Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=200, indent='hanging')
            for (cR in 1:length(listofcalls)) {
              tmpcall <- listofcalls[cR]
              Rmimic::typewriter(tmpcall, tabs=1, spaces=0, characters=200, indent='hanging')
            }
            
          } # something is being modeled
        } # DV select
      } # dataframe select
        
      invisible(shiny::stopApp())
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}
