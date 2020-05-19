#' pop_RmimicANOVA
#'
#' @description Popup window to compute SPSS style univariate ANOVA with effect size and confidence intervals
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, May 5, 2020
#'
#' @importFrom miniUI miniPage gadgetTitleBar miniContentPanel miniTitleBarCancelButton miniTitleBarButton
#' @importFrom shiny uiOutput renderUI wellPanel observeEvent stopApp runGadget dialogViewer icon textInput HTML
#' @importFrom shinyWidgets pickerInput actionBttn
#' @importFrom pkgcond suppress_conditions
#' 
#' @export

pop_RmimicANOVA <- function() {
  
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
    miniUI::gadgetTitleBar("Rmimic: Compute ANOVA", left = miniUI::miniTitleBarCancelButton(), right=NULL),
    
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
        plannedcontrasttext <- shiny::HTML("Use this field to indicate any planned contrasts.<br><small>Contrasts should be provided in quotes seperated by a comma, interactions should take the form variableX:variableY.</small>")
        
        # determine what is presented based upon format
        if (input$select_inputstyle == "Yes - Wide Format") {
          
          shiny::wellPanel(
            shiny::tagList(
              
              p('This function cannot presently handle wide format data. 
                Please transition your data into long format with your dependent
                variable all in a single column with independent variables in seperate columns. 
                There are a number of packages that can reshape a dataset from wide to long. 
                Another option is to subset dataframes and combine them vertically using the rbind function.')
              
            ) # end taglist
          ) # end wellPanel
          
        } else {
          
          shiny::wellPanel(
            shiny::tagList(
              shinyWidgets::pickerInput(
                inputId = "select_DV",
                label = "Select the variable of interest (dependent variable):",
                choices = workingdatavariables, width = "100%",
                multiple = FALSE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::pickerInput(
                inputId = "select_BSIV",
                label = "Select any between subject variables:",
                choices = workingdatavariables, width = "100%",
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
              ),
              
              shinyWidgets::pickerInput(
                inputId = "select_WSIV",
                label = "Select any within subject variables:",
                choices = workingdatavariables, width = "100%",
                multiple = TRUE,
                choicesOpt = list(subtext = workingdatavariablestypes),
                selected = NULL
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
                inputId = "select_sphericitystyle",
                label = "What type of sphericity correction should be run?",
                choices = c("Greenhouse-Geisser","Huynh-Feldt"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_teststyle",
                label = "What type of posthoc t-tests should be run?",
                choices = c("Parametric","Nonparametric"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "80%",
              ),
              
              shinyWidgets::radioGroupButtons(
                inputId = "select_posthocmethod",
                label = "What posthoc correction approach should be used?",
                choices = c("False Discovery Rate Control","Holm-Bonferroni", "Bonferroni", "Sidak", "Scheffe", "Tukey"),
                checkIcon = list(yes = shiny::icon("ok", lib = "glyphicon")),
                size = "normal",
                width = "90%",
              ),
              
              shiny::textInput("text_planned",
                       label = plannedcontrasttext,
                       value = "", 
                       width = "90%",
                       placeholder = NULL)
              
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
        
        # Check that selections are made
        if (!is.null(input$select_dataframe)) {
          
          if (input$select_inputstyle == "Yes - Wide Format") {
            # someday could implement this.
            
          } else {
            
            if (!is.null(input$select_DV)) {
              # user has chosen a DV
              if ((!is.null(input$select_BSIV)) | (!is.null(input$select_WSIV))) {
                # user has chosen at least either a between subjects or within subjects variable
                tmpcall <- 'result <- Rmimic::RmimicAnova('
                tmpcall <- sprintf('%sdata = %s', tmpcall, input$select_dataframe)
                tmpcall <- sprintf('%s, \n dependentvariable=%s', tmpcall, sprintf("'%s'",input$select_DV))
                if (!is.null(input$select_subID)) {
                  tmpcall <- sprintf('%s, \n subjectid=%s', tmpcall, sprintf("'%s'", input$select_subID[1]))
                } else {
                  tmpcall <- sprintf('%s, \n subjectid=NULL', tmpcall)
                }
                if (!is.null(input$select_BSIV)) {
                  tmpcall <- sprintf('%s, \n between=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_BSIV), collapse=", "))
                } else {
                  tmpcall <- sprintf('%s, \n between=NULL', tmpcall)
                }
                if (!is.null(input$select_WSIV)) {
                  tmpcall <- sprintf('%s, \n within=c(%s)', tmpcall, paste(sprintf("'%s'",input$select_WSIV), collapse=", "))
                } else {
                  tmpcall <- sprintf('%s, \n within=NULL', tmpcall)
                }
                if (input$select_sphericitystyle == 'Greenhouse-Geisser') {
                  tmpcall <- sprintf('%s, \n sphericity=%s', tmpcall, sprintf("'%s'", 'Greenhouse-Geisser'))
                } else {
                  tmpcall <- sprintf('%s, \n sphericity=%s', tmpcall, sprintf("'%s'", 'Huynh-Feldt'))
                }
                if (input$select_teststyle == 'Parametric') {
                  tmpcall <- sprintf('%s, \n nonparametric=FALSE', tmpcall)
                } else {
                  tmpcall <- sprintf('%s, \n nonparametric=TRUE', tmpcall)
                }
                if (input$select_posthocmethod == 'False Discovery Rate Control') {
                  tmpcall <- sprintf('%s, \n posthoc=%s', tmpcall, sprintf("'%s'", 'False Discovery Rate Control'))
                } else {
                  tmpcall <- sprintf('%s, \n posthoc=%s', tmpcall, sprintf("'%s'", input$select_posthocmethod))
                }
                if (input$text_planned != '') {
                  tmpcall <- sprintf('%s, \n planned=c(%s)', tmpcall, paste(sprintf("'%s'",input$text_planned), collapse=", "))
                } else {
                  tmpcall <- sprintf('%s, \n planned=NULL', tmpcall)
                }
                tmpcall <- sprintf('%s, \n feffect=%s, \n studywiseAlpha=0.05, confidenceinterval=0.95, verbose=TRUE)', tmpcall, sprintf("'%s'", "Generalized Eta Squared"))
                
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
                
              } # BS or WS variable chosen
            } # DV chosen
            
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
