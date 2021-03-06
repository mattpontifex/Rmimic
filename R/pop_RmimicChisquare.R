#' pop_RmimicChisquare
#'
#' @description Popup window to compute SPSS style regression analysis with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 7, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer HTML
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
  
  # Extract current tables from environment
  environelementstables <- ls(envir=.GlobalEnv) # global elements
  info_dfstables <- ""
  if (length(environelementstables) > 0) {
    environelementstablesidx <- which(sapply(environelementstables, function(x) is.table(get(x)))) # only tables
    environelementstables <- environelementstables[environelementstablesidx] 
    
    info_dfstables <- lapply(
      X = environelementstables,
      FUN = function(x) {
        tmp <- get(x, envir = .GlobalEnv)
        sprintf("table %d factors of  %d variables", nrow(tmp), ncol(tmp))
      }
    )
  }
  environelements <- c(environelements, environelementstables)
  info_dfs <- c(info_dfs, info_dfstables)
  
  
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
        if (!is.data.frame(workingdata)) {
          workingdata <- table2frame(workingdata)
        }
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
        
        model1contrasttext <- shiny::HTML("Only show significant breakdowns?<br><small>
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
            
            
            shinyWidgets::radioGroupButtons(
              inputId = "select_posthocmethod",
              label = "What posthoc correction approach should be used?",
              choices = c("False Discovery Rate Control","Holm-Bonferroni", "Bonferroni", "Sidak"),
              checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
              size = "normal",
              width = "90%",
            ),
            
            
          )
        )
      }
    })
    
    shiny::observeEvent(input$select_DV, {
      
      if (!is.null(input$select_dataframe)) {
        workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
        if (!is.data.frame(workingdata)) {
          workingdata <- table2frame(workingdata)
        }
        
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
        workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
        if (!is.data.frame(workingdata)) {
          workingdata <- table2frame(workingdata)
        }
        
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
    
    toListen <- shiny::reactive({
      list(input$done,input$generatecode)
    })
    
    shiny::observeEvent(toListen(), {
      if ((input$generatecode != 0) | (input$done != 0)) {
        
        # Check that selections are made
        if (!is.null(input$select_dataframe)) {
          
          if ((!is.null(input$select_DV)) & (!is.null(input$select_IV))) {
            
            # user has chosen a DV
            booltable <- FALSE
            workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
            if (!is.data.frame(workingdata)) {
              workingdata <- table2frame(workingdata)
              booltable <- TRUE
            }
            
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
            
            # array input
            #tmpcall <- sprintf('%svariables=c(%s, %s)', tmpcall, sprintf("'%s'",input$select_IV[1]), sprintf("'%s'",input$select_DV[1]))
            # specific input
            tmpcall <- sprintf('%sx=%s', tmpcall, sprintf("'%s'",input$select_IV[1]))
            tmpcall <- sprintf('%s, y=%s', tmpcall, sprintf("'%s'",input$select_DV[1]))
            
            tmpcall <- sprintf('%s, data=%s', tmpcall, input$select_dataframe)
            
            if (input$select_posthocmethod == 'False Discovery Rate Control') {
              tmpcall <- sprintf('%s, \n posthoc=%s', tmpcall, sprintf("'%s'", 'False Discovery Rate Control'))
            } else {
              tmpcall <- sprintf('%s, \n posthoc=%s', tmpcall, sprintf("'%s'", input$select_posthocmethod))
            }
            
            if (input$select_modelstyle == "Yes - Restrict output") {
              tmpcall <- sprintf('%s, planned=FALSE', tmpcall)
            } else {
              tmpcall <- sprintf('%s, planned=TRUE', tmpcall)
            }
            tmpcall <- sprintf('%s, \n studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE)', tmpcall)
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
              Rmimic::typewriter('Equivalent call:', tabs=0, spaces=0, characters=200, indent='hanging')
              codelevel <- 1
            }
            for (cR in 1:length(listofcalls)) {
              tmpcall <- listofcalls[cR]
              Rmimic::typewriter(tmpcall, tabs=codelevel, spaces=0, characters=200, indent='hanging')
            }
              
          } # DV and IV select
        } # dataframe select
        
        invisible(shiny::stopApp())
      }
    })
    shiny::observeEvent(input$cancel, {
      pkgcond::suppress_conditions(invisible(shiny::stopApp()))
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}