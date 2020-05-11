#' pop_imputation
#'
#' @description Popup window to impute missing data
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 10, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#'
#' @export

pop_imputation <- function() {
  
  # Get the document context.
  #context <- rstudioapi::getActiveDocumentContext()
  
  # Extract current data from environment
  environelements <- ls(envir=.GlobalEnv) # global elements
  environelementsidx <- which(sapply(environelements, function(x) is.data.frame(get(x)))) # only dataframes
  environelements <- environelements[environelementsidx] 
  computetext <- HTML("Compute<br><small>
                                   Note: This may take a few moments.</small>")
  
  info_dfs <- lapply(
    X = environelements,
    FUN = function(x) {
      tmp <- get(x, envir = .GlobalEnv)
      sprintf("%d obs. of  %d variables", nrow(tmp), ncol(tmp))
    }
  )
  info_dfs <- unlist(info_dfs)
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Rmimic: Impute Missing Observations", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
      label = computetext,
      style = "simple", 
      color = "primary",
      block=TRUE,
      size="md"
    )
    
  )
  
  server <- function(input, output, session) {
    
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
          if (sum(!complete.cases(testvect)) > 0) {
            workingdatavariablestypes <- c(workingdatavariablestypes, sprintf('    type: %s,       observations missing: %d', testvectout, sum(!complete.cases(testvect))))
          } else {
            workingdatavariablestypes <- c(workingdatavariablestypes, sprintf('    type: %s', testvectout))
          }
        }
        
        methodoptions <- c("Predictive mean matching",
                           "Weighted predictive mean matching", 
                           "Random sample from observed values",
                           "Unconditional mean imputation",
                           "Logistic regression")
        imputationtext <- HTML("Select the method for imputation:<br><small>
                                   See mice help file for a complete list of methods.</small>")
        imputationnumtext <- HTML("How many imputations should be performed?<br><small>
                                   The larger the number of imputations, the longer this will take.</small>")
        
        
        shiny::wellPanel(
          shiny::tagList(
            shinyWidgets::pickerInput(
              inputId = "select_variables",
              label = "Select variables to impute:",
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            br(),
            
            shinyWidgets::pickerInput(
              inputId = "select_method",
              label = imputationtext,
              choices = methodoptions, width = "80%",
              options = list(), 
              multiple = FALSE,
              selected = 1
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_restriction",
              label = "Should imputed values be restricted to range of observed values?",
              choices = c("Yes","No"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "80%",
            ),
            
            shinyWidgets::sliderTextInput(
              inputId = "select_numofimputations",
              label = imputationnumtext, 
              choices = seq(from = 2, to = 50, by = 2),
              selected = 10,
              width = "80%"
            )
          )
        )
      }
    })
    
    shiny::observeEvent(input$done, {
        
      # Check that selections are made
      if (!is.null(input$select_dataframe)) {
        
        tmpcall <- sprintf('%s <- Rmimic::multipleimputation(data=%s',input$select_dataframe,input$select_dataframe)
        
        if (!is.null(input$select_variables)) {
          workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
          # Some variables were chosen
          workingdatavariables <- names(workingdata)
          if (paste(sprintf("'%s'",input$select_variables), collapse=", ") == paste(sprintf("'%s'",workingdatavariables), collapse=", ")) {
            # all variables were selected
            
          } else {
            # specific variables were selected
            tmpcall <- sprintf('%s, variables=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_variables), collapse=", "))
          }
        }
        
        if (input$select_method == "Predictive mean matching") {
          tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'", 'pmm'))
        } else if (input$select_method == "Weighted predictive mean matching") {
          tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'", 'midastouch'))
        } else if (input$select_method == "Random sample from observed values") {
          tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'", 'sample'))
        } else if (input$select_method == "Unconditional mean imputation") {
          tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'", 'mean'))
        } else if (input$select_method == "Logistic regression") {
          tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'", 'logreg'))
        }
        
        if (input$select_restriction == "Yes") {
          tmpcall <- sprintf('%s, restrict=TRUE', tmpcall)
        } else {
          tmpcall <- sprintf('%s, restrict=FALSE', tmpcall)
        }
        tmpcall <- sprintf('%s, imputations=%d)', tmpcall, input$select_numofimputations)
        
        # execute call
        boolattempt <- FALSE
        boolattempt <- tryCatch({
          eval(parse(text=tmpcall))
          boolattempt <- TRUE}
        )
        if (boolattempt == FALSE) {
          Rmimic::typewriter('Uh oh.. Something went wrong. But the syntax for the function is provided below.', tabs=0, spaces=0, characters=80, indent='hanging')
        }
        
        Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=200, indent='hanging')
        Rmimic::typewriter(tmpcall, tabs=1, spaces=0, characters=200, indent='hanging')
      
      }
      
      invisible(shiny::stopApp())
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}

