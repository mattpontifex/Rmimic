#' pop_identifyoutliers
#'
#' @description Popup window to impute missing data
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 18, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer HTML
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#' @export

pop_identifyoutliers <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Identify and Remove Outliers", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
        
        iqrlimitnumtext <- shiny::HTML("Select the interquartile limit:<br><small>
                                   SPSS identifies an interquartile limit of 3 times the interquartile range as an extreme outlier.</small>")
        
        shiny::wellPanel(
          shiny::tagList(
            shinyWidgets::pickerInput(
              inputId = "select_variables",
              label = "Select variables to identify outliers within:",
              choices = workingdatavariables, width = "100%",
              options = list(`actions-box` = TRUE),
              multiple = TRUE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shiny::br(),
            
            shinyWidgets::sliderTextInput(
              inputId = "select_iqrlimit",
              label = iqrlimitnumtext, 
              choices = seq(from = 1, to = 5, by = 0.25),
              selected = 3,
              width = "80%"
            ),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_restriction",
              label = "Should identification of outliers be restricted to a particular tail?",
              choices = c("Two tailed", "Upper only", "Lower only"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "80%"
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
          
          tmpcall <- sprintf('%s <- Rmimic::identifyoutliersdataframe(data=%s',input$select_dataframe,input$select_dataframe)
          
          if (!is.null(input$select_variables)) {
            workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
            # Some variables were chosen
            workingdatavariables <- names(workingdata)
            if (paste(sprintf("'%s'",input$select_variables), collapse=", ") == paste(sprintf("'%s'",workingdatavariables), collapse=", ")) {
              # all variables were selected
              
            } else {
              # specific variables were selected
              tmpcall <- sprintf('%s, variables=c(%s) \n', tmpcall, paste(sprintf("'%s'",input$select_variables), collapse=", "))
            }
          }
          tmpcall <- sprintf('%s, iqrlimit=%d', tmpcall, input$select_iqrlimit)
          
          
          if (input$select_restriction == "Two tailed") {
            tmpcall <- sprintf('%s, direction=%s', tmpcall, sprintf("'%s'", 'both'))
          } else if (input$select_restriction == "Upper only") {
            tmpcall <- sprintf('%s, direction=%s', tmpcall, sprintf("'%s'", 'upperonly'))
          } else if (input$select_restriction == "Lower only") {
            tmpcall <- sprintf('%s, direction=%s', tmpcall, sprintf("'%s'", 'loweronly'))
          }
          tmpcall <- sprintf('%s, verbose=TRUE)', tmpcall)
          
          # execute call
          codelevel <- 0 
          if (input$done) {
            boolattempt <- FALSE
            boolattempt <- tryCatch({
              eval(parse(text=tmpcall), envir = .GlobalEnv)
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
        
        invisible(shiny::stopApp())
      }
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}
