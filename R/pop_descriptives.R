#' pop_descriptives
#'
#' @description Popup window to compute SPSS style descriptives.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 3, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#'
#' @export

pop_descriptives <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Compute Descriptives", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
        shiny::wellPanel(
        shiny::tagList(
          shinyWidgets::pickerInput(
            inputId = "select_variables",
            label = "Select/deselect the variables of interest:",
            choices = workingdatavariables, width = "100%",
            options = list(`actions-box` = TRUE),
            multiple = TRUE,
            choicesOpt = list(subtext = workingdatavariablestypes),
            selected = NULL
          ),
          
          shinyWidgets::pickerInput(
            inputId = "select_groupvariables",
            label = "Seperate by:",
            choices = workingdatavariables, width = "100%",
            options = list(), 
            multiple = TRUE,
            choicesOpt = list(subtext = workingdatavariablestypes),
            selected = NULL
          )
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
          if (!is.null(input$select_variables)) {
            tmpcall <- 'desc <- Rmimic::descriptives('
            # Extract current data from environment
            workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
            tmpcall <- sprintf('%svariables=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_variables), collapse=", "))
            if (!is.null(input$select_groupvariables)) {
              tmpcall <- sprintf('%s, groupvariable=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_groupvariables), collapse=", "))
            }
            tmpcall <- sprintf('%s, data=%s', tmpcall, input$select_dataframe)
            tmpcall <- sprintf('%s, verbosedescriptives=TRUE, verbosefrequencies=TRUE)', tmpcall)
            
            # execute call
            codelevel <- 0 
            if (input$done) {
              boolattempt <- FALSE
              boolattempt <- tryCatch({
                eval(parse(text=tmpcall))
                boolattempt <- TRUE}
              )
              if (boolattempt == FALSE) {
                Rmimic::typewriter('Uh oh.. Something went wrong. But the syntax for the function is provided below.', tabs=0, spaces=0, characters=80, indent='hanging')
              }
              Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=200, indent='hanging')
              codelevel <- 1
            }
            Rmimic::typewriter(tmpcall, tabs=codelevel, spaces=0, characters=200, indent='hanging')
          }
        }
        
        invisible(shiny::stopApp())
      }
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}
