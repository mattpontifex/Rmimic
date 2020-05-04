#' pop_descriptives
#'
#' @description Popup window to compute SPSS style descriptives.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 3, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny selectInput uiOutput renderUI wellPanel observeEvent stopApp runGadget
#' 
#'
#' @export

pop_descriptives <- function() {
  
  # Get the document context.
  #context <- rstudioapi::getActiveDocumentContext()
  
  # Extract current data from environment
  environelements <- ls(envir=.GlobalEnv) # global elements
  environelementsidx <- which(sapply(environelements, function(x) is.data.frame(get(x)))) # only dataframes
  environelements <- environelements[environelementsidx] 
  environelementsidx <- unlist(as.numeric(1:length(environelementsidx)))
  names(environelementsidx) <- environelements
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Rmimic: Compute Descriptives", left = miniUI::miniTitleBarCancelButton(),
                   right = miniUI::miniTitleBarButton("done", "Compute", primary = TRUE)),
    miniUI::miniContentPanel(
      shiny::selectInput("select_dataframe", label = "Select your working data frame:", 
              choices = environelementsidx, 
              selected = 1),
      shiny::uiOutput("ui1"), # This outputs the dynamic UI component
    )
  )
  
  server <- function(input, output) {
    
    output$ui1 <- shiny::renderUI({
      if (is.null(input$select_dataframe)) {
        return()
      } else {
        
        # Extract current data from environment
        environelements <- ls(envir=.GlobalEnv) # global elements
        environelementsidx <- which(sapply(environelements, function(x) is.data.frame(get(x)))) # only dataframes
        environelements <- environelements[environelementsidx] 
        environelementsidx <- unlist(as.numeric(1:length(environelementsidx)))
        
        workingdata <- get(environelements[as.numeric(input$select_dataframe)], envir = .GlobalEnv)
        workingdatavariables <- names(workingdata)
        names(workingdatavariables) <- workingdatavariables
        shiny::wellPanel(
          shiny::selectInput("select_variables", "Select the variables:",
            choices = workingdatavariables, selected = NULL, multiple = TRUE
          ),
          shiny::selectInput("select_groupvariables", "Seperate by:",
                      choices = workingdatavariables, selected = NULL, multiple = TRUE
          )
        )
      }
    })
    
    shiny::observeEvent(input$done, {
      
      # Check that selections are made
      if (!is.null(input$select_dataframe)) {
        if (!is.null(input$select_variables)) {
          tmpcall <- 'desc <- descriptives('
          
          # Extract current data from environment
          environelements <- ls(envir=.GlobalEnv) # global elements
          environelementsidx <- which(sapply(environelements, function(x) is.data.frame(get(x)))) # only dataframes
          environelements <- environelements[environelementsidx] 
          environelementsidx <- unlist(as.numeric(1:length(environelementsidx)))
          workingdataname <- environelements[as.numeric(input$select_dataframe)]
          workingdata <- get(workingdataname, envir = .GlobalEnv)
          tmpcall <- sprintf('%svariables=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_variables), collapse=", "))
          if (!is.null(input$select_groupvariables)) {
            tmpcall <- sprintf('%s, groupvariable=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_groupvariables), collapse=", "))
          }
          tmpcall <- sprintf('%s, data=%s', tmpcall, workingdataname)
          tmpcall <- sprintf('%s, verbosedescriptives=TRUE, verbosefrequencies=TRUE)', tmpcall)
          
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
  
  shiny::runGadget(ui, server, viewer = dialogViewer(""))
}