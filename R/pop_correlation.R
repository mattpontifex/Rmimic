#' pop_correlation
#'
#' @description Popup window to compute SPSS style correlations
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 4, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer
#' @importFrom shinyWidgets pickerInput actionBttn
#' 
#'
#' @export

pop_correlation <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Compute Correlations", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
        
        methodoptions <- c("Pearson", "Spearman", "Kendall")
        methodoptionsnotes <- c("product moment correlation coefficient", "rho correlation coefficients", "tau correlation coefficients")
        
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
            inputId = "select_partials",
            label = "Partial by:",
            choices = workingdatavariables, width = "100%",
            options = list(), 
            multiple = TRUE,
            choicesOpt = list(subtext = workingdatavariablestypes),
            selected = NULL
          ),
          
          
          shinyWidgets::pickerInput(
            inputId = "select_method",
            label = "Select correlation coefficient:",
            choices = methodoptions, width = "80%",
            options = list(), 
            multiple = FALSE,
            choicesOpt = list(subtext = methodoptionsnotes),
            selected = 1
          ),
          
        )
          
        )
      }
    })
    
    shiny::observeEvent(input$done, {
      
      # Check that selections are made
      if (!is.null(input$select_dataframe)) {
        if (!is.null(input$select_variables)) {
          tmpcall <- 'sR <- Rmimic::correlation('
          # Extract current data from environment
          workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
          tmpcall <- sprintf('%svariables=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_variables), collapse=", "))
          if (!is.null(input$select_partials)) {
            tmpcall <- sprintf('%s, partial=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_partials[1]), collapse=", "))
          }
          tmpcall <- sprintf('%s, data=%s', tmpcall, input$select_dataframe)
          if (!is.null(input$select_method)) {
            tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'%s'",tolower(input$select_method)))
          } else {
            tmpcall <- sprintf('%s, method=%s', tmpcall, sprintf("'pearson'"))
          }
          tmpcall <- sprintf('%s, studywiseAlpha=0.05, confidenceinterval=0.95, listwise=TRUE, verbose=TRUE)', tmpcall)
          
          # execute call
          eval(parse(text=tmpcall))
          Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=80, indent='hanging')
          Rmimic::typewriter(tmpcall, tabs=1, spaces=0, characters=80, indent='hanging')
        }
      }
      
      invisible(shiny::stopApp())
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer(""))
}
