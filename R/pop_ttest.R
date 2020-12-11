#' pop_ttest
#'
#' @description Popup window to compute SPSS style t-test with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 4, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer icon
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#' @export

pop_ttest <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Compute T-test", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
      
      shinyWidgets::radioGroupButtons(
        inputId = "select_inputstyle",
        label = "Is your dependent variable in seperate columns?",
        choices = c("No - Long Format","Yes - Wide Format"),
        checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
        size = "normal"
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
        
        # determine what is presented based upon format
        
        if (input$select_inputstyle == "Yes - Wide Format") {
          
          shiny::wellPanel(
            shiny::tagList(
              shinyWidgets::pickerInput(
                inputId = "select_variablesWF",
                label = "Select the variables for comparison:",
                choices = workingdatavariables, width = "100%",
                options = list(`actions-box` = TRUE),
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_IVstyleWF",
                label = "What type of comparison is the data?",
                choices = c("Between Subjects","Within Subjects"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal"
              ),
              
              shinyWidgets::pickerInput(
                inputId = "select_subIDWF",
                label = "Which variable corresponds to the subject ID (if available)?",
                choices = workingdatavariables, width = "100%",
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_teststyleWF",
                label = "What type of tests should be run?",
                choices = c("Parametric","Nonparametric"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_posthocmethodWF",
                label = "What posthoc tests should be run?",
                choices = c("Holm-Bonferroni", "Bonferroni", "Sidak"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
            ) # end taglist
          ) # end wellPanel
          
        } else {
          
          shiny::wellPanel(
            shiny::tagList(
              shinyWidgets::pickerInput(
                inputId = "select_DV",
                label = "Select the variables of interest (dependent variables):",
                choices = workingdatavariables, width = "100%",
                options = list(`actions-box` = TRUE),
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::pickerInput(
                inputId = "select_IV",
                label = "Select the independent variable:",
                choices = workingdatavariables, width = "100%",
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_IVstyle",
                label = "Is your independent variable:",
                choices = c("Between Subjects","Within Subjects"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal"
              ),
              
              shinyWidgets::pickerInput(
                inputId = "select_subID",
                label = "Which variable corresponds to the subject ID (if available)?",
                choices = workingdatavariables, width = "100%",
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_teststyle",
                label = "What type of tests should be run?",
                choices = c("Parametric","Nonparametric"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_posthocmethod",
                label = "What posthoc tests should be run?",
                choices = c("Holm-Bonferroni", "Bonferroni", "Sidak"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
            )
          )
        } # end input style
        
      }
    })
    
    toListen <- shiny::reactive({
      list(input$done,input$generatecode)
    })
    
    shiny::observeEvent(toListen(), {
      if ((input$generatecode != 0) | (input$done != 0)) {
        
        listofcalls <- c()
        
        # Check that selections are made
        if (!is.null(input$select_dataframe)) {
          
          if (input$select_inputstyle == "Yes - Wide Format") {
          
            if (!is.null(input$select_variablesWF)) {
              
              # show user how to take data from wide to long format
              tmpcall <- sprintf('workingdata <- rbind(')
              for (cR in 1:length(input$select_variablesWF)) {
                tmpcall <- sprintf('%sdata.frame(%s=%s', tmpcall, sprintf("'%s'", 'DV'), input$select_dataframe)
                tmpcall <- sprintf('%s$%s',tmpcall, input$select_variablesWF[cR])
                if (!is.null(input$select_subIDWF)) {
                  tmpcall <- sprintf('%s,%s=%s$%s', tmpcall, sprintf("'%s'", input$select_subIDWF[1]), input$select_dataframe, input$select_subIDWF[1])
                }
                tmpcall <- sprintf('%s,%s=%s)', tmpcall, sprintf("'%s'", 'IV'), sprintf("'%s'", input$select_variablesWF[cR]))
                if (cR < (length(input$select_variablesWF))) {
                  tmpcall <- sprintf('%s, ',tmpcall)
                }
              }
              tmpcall <- sprintf('%s)',tmpcall)
              tmpcallseg1 <- tmpcall
              #eval(parse(text=tmpcall))
              listofcalls <- c(listofcalls, tmpcall)
              
              tmpcall <- sprintf('ttestresult <- Rmimic::RmimicTtest(workingdata, dependentvariable=%s', sprintf("'DV'"))
               if (!is.null(input$select_subIDWF)) {
                tmpcall <- sprintf('%s, \n subjectid=%s', tmpcall, sprintf("'%s'", input$select_subIDWF[1]))
              } else {
                tmpcall <- sprintf('%s, subjectid=NULL', tmpcall)
              }
              if (input$select_IVstyleWF == 'Between Subjects') {
                tmpcall <- sprintf('%s, \n between=c(%s)', tmpcall, sprintf("'%s'", 'IV'))
                tmpcall <- sprintf('%s, within=NULL', tmpcall)
              } else {
                tmpcall <- sprintf('%s, between=NULL', tmpcall)
                tmpcall <- sprintf('%s, \n within=c(%s)', tmpcall, sprintf("'%s'", 'IV'))
              }
              if (input$select_teststyleWF == 'Parametric') {
                tmpcall <- sprintf('%s, \n nonparametric=FALSE', tmpcall)
              } else {
                tmpcall <- sprintf('%s, \n nonparametric=TRUE', tmpcall)
              }
              if (input$select_posthocmethodWF == 'None') {
                tmpcall <- sprintf('%s, posthoc=FALSE', tmpcall)
              } else {
                tmpcall <- sprintf('%s, posthoc=%s', tmpcall, sprintf("'%s'", input$select_posthocmethodWF))
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
            }
            
          } else {
            if ((!is.null(input$select_DV)) & (!is.null(input$select_IV))) {
              tmpcall <- 'ttestresult <- Rmimic::RmimicTtest('
              # Extract current data from environment
              workingdata <- get(input$select_dataframe, envir = .GlobalEnv)
              tmpcall <- sprintf('%s%s', tmpcall, input$select_dataframe)
              tmpcall <- sprintf('%s, \n dependentvariable=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_DV), collapse=", "))
              if (!is.null(input$select_subID)) {
                tmpcall <- sprintf('%s, \n subjectid=%s', tmpcall, sprintf("'%s'", input$select_subID[1]))
              } else {
                tmpcall <- sprintf('%s, subjectid=NULL', tmpcall)
              }
              if (input$select_IVstyle == 'Between Subjects') {
                tmpcall <- sprintf('%s, \n between=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_IV), collapse=", "))
                tmpcall <- sprintf('%s, within=NULL', tmpcall)
              } else {
                tmpcall <- sprintf('%s, between=NULL', tmpcall)
                tmpcall <- sprintf('%s, \n within=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_IV), collapse=", "))
              }
              if (input$select_teststyle == 'Parametric') {
                tmpcall <- sprintf('%s, \n nonparametric=FALSE', tmpcall)
              } else {
                tmpcall <- sprintf('%s, \n nonparametric=TRUE', tmpcall)
              }
              if (input$select_posthocmethod == 'None') {
                tmpcall <- sprintf('%s, posthoc=FALSE', tmpcall)
              } else {
                tmpcall <- sprintf('%s, posthoc=%s', tmpcall, sprintf("'%s'", input$select_posthocmethod))
              }
              tmpcall <- sprintf('%s, \n studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE)', tmpcall)
              
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
            
          } # end inputstyle
        } #end empty dataframe
        
        invisible(shiny::stopApp())
      }
    })
    shiny::observeEvent(input$cancel, {
      invisible(shiny::stopApp())
    })
  }
  
  pkgcond::suppress_conditions(shiny::runGadget(ui, server, viewer = shiny::dialogViewer("")))
}
