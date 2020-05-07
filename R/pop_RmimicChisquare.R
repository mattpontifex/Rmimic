#' pop_RmimicChisquare
#'
#' @description Popup window to compute SPSS style regression analysis with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 7, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#'
#' @export

pop_RmimicChisquare <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Compute Chi-square", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
          workingdatavariablestypes <- c(workingdatavariablestypes, sprintf('    type: %s', testvectout))
        }
        
        model1contrasttext <- HTML("Only show significant breakdowns?<br><small>
                                   Results of the overall model will still be shown.</small>")
        
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
            
            
            shiny::textInput("text_factorDV",
                             label = 'Specify the factor levels of the Outcome variable:',
                             value = "", 
                             width = "80%",
                             placeholder = NULL),
            
            br(),
            br(),
            
            shinyWidgets::pickerInput(
              inputId = "select_IV",
              label = 'Select the Predictor variable (independent variable)',
              choices = workingdatavariables, width = "100%",
              multiple = FALSE,
              choicesOpt = list(subtext = workingdatavariablestypes),
              selected = NULL
            ),
            
            shiny::textInput("text_factorIV",
                             label = 'Specify the factor levels of the Predictor variable:',
                             value = "", 
                             width = "80%",
                             placeholder = NULL),
            
            br(),
            br(),
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_modelstyle",
              label = model1contrasttext,
              choices = c("Yes - Restrict output","No - Show me everything"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "80%",
            ),
            
          )
        )
      }
    })
    
    shiny::observeEvent(input$select_DV, {
      
      if (!is.null(input$select_dataframe)) {
        if (!is.null(input$select_DV)) {
          name <- levels(workingdata[,input$select_DV[1]])
          if (is.null(name[1])) {
            name <- unique(as.character(workingdata[,input$select_DV[1]]))
          }
          
          name <- paste(sprintf("'%s'", as.character(name)), sep="' '", collapse=", ")
          shiny::updateTextInput(session, "text_factorDV", value=name)
        }
      }
    })
    
    
    shiny::observeEvent(input$select_IV, {
      
      if (!is.null(input$select_dataframe)) {
        if (!is.null(input$select_IV)) {
          name <- levels(workingdata[,input$select_IV[1]])
          if (is.null(name[1])) {
            name <- unique(as.character(workingdata[,input$select_IV[1]]))
          }
          
          name <- paste(sprintf("'%s'", as.character(name)), sep="' '", collapse=", ")
          shiny::updateTextInput(session, "text_factorIV", value=name)
        }
      }
    })
    
    shiny::observeEvent(input$done, {
      
      # Check that selections are made
      if (!is.null(input$select_dataframe)) {
        if ((!is.null(input$select_DV)) & (!is.null(input$select_IV))) {
          # user has chosen a DV
            
          listofcalls <- c()
          
          boolrep <- FALSE
          name <- levels(workingdata[,input$select_DV[1]])
          if (is.null(name[1])) {
            boolrep <- TRUE # needs to be factored
            name <- unique(as.character(workingdata[,input$select_DV[1]]))
          }
          name <- paste(sprintf("'%s'", as.character(name)), sep="' '", collapse=", ")
          comptext <- as.character(input$text_factorDV)
          if (name != comptext) {
            boolrep <- TRUE # needs to be factored
            name <- comptext
          }
          if (boolrep == TRUE) {
            tmpcall <- sprintf('%s[,%s] <- factor(%s[,%s]', input$select_dataframe, sprintf("'%s'",input$select_DV[1]), input$select_dataframe, sprintf("'%s'",input$select_DV[1]))
            tmpcall <- sprintf('%s, levels=c(%s))', tmpcall, name)
            listofcalls <- c(listofcalls, tmpcall)
          }
          
          boolrep <- FALSE
          name <- levels(workingdata[,input$select_IV[1]])
          if (is.null(name[1])) {
            boolrep <- TRUE # needs to be factored
            name <- unique(as.character(workingdata[,input$select_IV[1]]))
          }
          name <- paste(sprintf("'%s'", as.character(name)), sep="' '", collapse=", ")
          comptext <- as.character(input$text_factorIV)
          if (name != comptext) {
            boolrep <- TRUE # needs to be factored
            name <- comptext
          }
          if (boolrep == TRUE) {
            tmpcall <- sprintf('%s[,%s] <- factor(%s[,%s]', input$select_dataframe, sprintf("'%s'",input$select_IV[1]), input$select_dataframe, sprintf("'%s'",input$select_IV[1]))
            tmpcall <- sprintf('%s, levels=c(%s))', tmpcall, name)
            listofcalls <- c(listofcalls, tmpcall)
          }
          
          tmpcall <- 'chisquareresult <- Rmimic::RmimicChisquare('
          tmpcall <- sprintf('%svariables=c(%s, %s)', tmpcall, sprintf("'%s'",input$select_IV[1]), sprintf("'%s'",input$select_DV[1]))
          tmpcall <- sprintf('%s, data=%s', tmpcall, input$select_dataframe)
          if (input$select_modelstyle == "Yes - Restrict output") {
            tmpcall <- sprintf('%s, planned=FALSE', tmpcall)
          } else {
            tmpcall <- sprintf('%s, planned=TRUE', tmpcall)
          }
          tmpcall <- sprintf('%s, studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE)', tmpcall)
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
            
        } # DV and IV select
      } # dataframe select
      
      invisible(shiny::stopApp())
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}