#' pop_RmimicLMcontrast
#'
#' @description Popup window to compute SPSS style regression analysis with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 5, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer HTML reactive
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' @importFrom stats glm lm binomial
#' 
#'
#' @export

pop_RmimicLMcontrast <- function() {
  
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
      inputId = "generatecode",
      label = "Generate Code Only",
      style = "simple", 
      color = "royal",
      block=TRUE,
      size="sm"
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
        workingdata <- Rmimic::antitibbler(workingdata)
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
        
        model1contrasttext <- shiny::HTML("Model 1: Select the explanatory variables (independent variables) for the base model:<br><small>
                                   This can be a constant or control variables.</small>")
        model2contrasttext <- shiny::HTML("Model 2: Select the explanatory variables (independent variables) that should be added to the base model:")
        
        shiny::wellPanel(
          shiny::tagList(
            
            shinyWidgets::pickerInput(
              inputId = "select_DV",
              label = "Select the variable you are trying to predict (dependent variable):",
              choices = workingdatavariables, width = "100%",
              multiple = FALSE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = 0
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
              choices = c("All Variables Entered","Forward Stepwise","Bidirectional Stepwise"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "100%",
            ),
            
            
            shinyWidgets::pickerInput(
              inputId = "select_M2variables",
              label = model2contrasttext,
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_M2style",
              label = "How should Model 2 be run?",
              choices = c("All Variables Entered","Forward Stepwise","Bidirectional Stepwise"),
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
            # user has chosen a DV
            if ((!is.null(input$select_M1variables)) | (!is.null(input$select_M2variables))) {
              
              # basic housekeeping
              M1variables <- input$select_M1variables
              if (length(which(M1variables == '< constant >')) > 0) {
                M1variables[which(M1variables == '< constant >')] <- '1'
              }
              M2variables <- input$select_M2variables
              if (length(which(M2variables == '< constant >')) > 0) {
                M2variables[which(M2variables == '< constant >')] <- '1'
              }
              
              listofcalls <- c()
              tmpcall <- sprintf('%s <- Rmimic::antitibbler(%s)\n', input$select_dataframe, input$select_dataframe)
              listofcalls <- c(listofcalls, tmpcall)
              tmppref <- ''
              tempsuff <- ''
              if (input$select_modelstyle == 'Linear Regression') {
                tmppref <- 'basefit <- stats::lm('
                tmppref2 <- 'fit <- stats::lm('
                tempsuff <- ''
              } else {
                tmppref <- 'basefit <- stats::glm('
                tmppref2 <- 'fit <- stats::glm('
                tempsuff <- 'family=stats::binomial(link = "logit"), '
              }
              tempsuff <- sprintf('%sdata=%s)', tempsuff, input$select_dataframe)
              
              # Figure out what the user wants
              if ((length(M1variables) > 0) & (length(M2variables) > 0)) {
                # user entered both set of variables
                
                tmpformsimple <- sprintf('%s ~ 1', input$select_DV[1])
                
                # Model 1
                # populate model
                tmpform <- sprintf('%s ~ %s', input$select_DV[1], paste(M1variables, collapse=" + ")) 
                
                # select model approach
                modeloptions <- input$select_M1style
                if ((input$select_M1style == "Forward Stepwise") | (input$select_M1style == "Bidirectional Stepwise")) {
                  if ((length(which(M1variables == '1')) > 0) & (length(M1variables) == 1)) {
                    # only a constant was selected so cannot do stepwise
                    modeloptions <- "All Variables Entered"
                  }
                }
                
                # select model approach
                if (modeloptions == "All Variables Entered") {
                  tmpcall <- sprintf('%s%s, %s',tmppref,tmpform,tempsuff) 
                  listofcalls <- c(listofcalls, tmpcall)
                } else {
                  # start with constant
                  tmpcall <- sprintf('%s%s, %s',tmppref,tmpformsimple,tempsuff) 
                  listofcalls <- c(listofcalls, tmpcall)
                  if (modeloptions == "Forward Stepwise") {
                    tmpcall <- sprintf('basefit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="forward", scope=list(lower=%s, \n upper=%s), trace=TRUE)$terms)', tmpformsimple, tmpform)
                  } else {
                    tmpcall <- sprintf('basefit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="both", scope=list(lower=%s, \n upper=%s), trace=TRUE)$terms)', tmpformsimple, tmpform)
                  } 
                  listofcalls <- c(listofcalls, tmpcall)
                }
                
                # Model 2 - additive with model 1
                M2variables <- unique(c(M1variables, M2variables))
                M2variables <- M2variables[which(M2variables != '1')]
                
                # populate model
                tmpform <- sprintf('%s ~ %s', input$select_DV[1], paste(M2variables, collapse=" + ")) 
                
                # select model approach
                if (input$select_M2style == "All Variables Entered") {
                  tmpcall <- sprintf('%s%s, %s',tmppref2,tmpform,tempsuff) 
                  listofcalls <- c(listofcalls, tmpcall)
                } else {
                  if (input$select_M2style == "Forward Stepwise") {
                    tmpcall <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="forward", scope=list(upper=%s), trace=TRUE)$terms)', tmpform)
                  } else {
                    tmpcall <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="both", scope=list(lower=%s, \n upper=%s), trace=TRUE)$terms)', tmpformsimple, tmpform)
                  } 
                  listofcalls <- c(listofcalls, tmpcall)
                }
              
              } else {
                # user only entered one set of variables
                
                # create constant model
                tmpformsimple <- sprintf('%s ~ 1', input$select_DV[1])
                tmpcall <- sprintf('%s%s, %s',tmppref,tmpformsimple,tempsuff) 
                listofcalls <- c(listofcalls, tmpcall)
                
                # create second model
                if (length(M1variables) > 0) {
                  workingvariables <- M1variables
                  modeloptions <- input$select_M1style
                } else {
                  workingvariables <- M2variables
                  modeloptions <- input$select_M2style
                }
                tmpform <- sprintf('%s ~ %s', input$select_DV[1], paste(workingvariables, collapse=" + ")) 
                
                if (modeloptions == "All Variables Entered") {
                  tmpcall <- sprintf('fit <- stats::update(basefit,formula=%s)', tmpform) # update with new call
                } else {
                  if (modeloptions == "Forward Stepwise") {
                    tmpcall <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="forward", scope=list(lower=%s, \n upper=%s), trace=TRUE)$terms)', tmpformsimple, tmpform)
                  } else {
                    tmpcall <- sprintf('fit <- stats::update(basefit,formula=MASS::stepAIC(basefit, \n direction="both", scope=list(lower=%s, \n upper=%s), trace=TRUE)$terms)', tmpformsimple, tmpform)
                  } 
                }
                listofcalls <- c(listofcalls, tmpcall)
              }
              
              tmpcall <- 'regresult <- Rmimic::RmimicLMcontrast(basefit, fit, studywiseAlpha=0.05, confidenceinterval=0.95)'
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
              
            } # something is being modeled
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
