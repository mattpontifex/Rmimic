#' typewriter
#'
#' @description Function to control the text outputted to the console to create a consistent indent.
#'
#' @usage
#' typewriter(textstring, tabs, spaces, characters, indent)
#'
#' @param textstring String to print to console.
#' @param tabs How many tabs should be printed before the text.
#' @param spaces Following tab indents, how many spaces should be printed before the text.
#' @param characters How many characters wide should the output be. Default is options()$width .
#' @param indent Parameter to control the indent. Positive numbers indicate an indent. Negative numbers indicate hanging format.
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, August 31, 2017
#'
#'
#' @examples
#'     # Example text block
#'     textblock <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#'                   sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#'                   Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris
#'                   nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in
#'                   reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla
#'                   pariatur. Excepteur sint occaecat cupidatat non proident, sunt in
#'                   culpa qui officia deserunt mollit anim id est laborum."
#'
#'     # Output indented text 3 tabs in
#'     typewriter(textblock, tabs=3, spaces=0, characters=100, indent=TRUE)
#'
#'     # Output hanging text 1 tabs in with 70 character width
#'     typewriter(textblock, tabs=1, spaces=0, characters=70, indent='hanging')
#'
#'     # Output indented text 2 tabs and 3 spaces in, with 100 character width
#'     typewriter(textblock, tabs=2, spaces=3, characters=100, indent=FALSE)
#'
#'
#' @export

typewriter <- function(textstring = NULL, tabs = NULL, spaces = NULL, characters = NULL, indent=NULL) {
# revised 10-1-2017 - added ability to handle text that does not wrap and text with no spaces.

  wid <- floor(options()$width *.98)
  if (wid > 93) {
    wid <- 93
  }
  
  if (!is.null(textstring)) {

    if (is.null(indent)) {
      indent <- 0
    } else if (toupper(indent) == toupper("indent")) {
      indent <- 2
    } else if (toupper(indent) == toupper("hanging")) {
      indent <- -4
    } else if (indent == TRUE) {
      indent <- 2
    }  else {
      indent <- as.numeric(indent)
    }

    if (is.null(tabs)) {
      tabs <- 0
    } else {
      tabs <- as.numeric(tabs)
    }
    if (is.null(spaces)) {
      spaces <- 0
    } else {
      spaces <- as.numeric(spaces)
    }

    if (is.null(characters)) {
      characters <- as.numeric(options()$width)
    } else {
      characters <- as.numeric(characters)
    }

    # Create the first line spacing
    firstlinespacing <- ""
    if ((!is.null(tabs)) & (tabs != 0)) {
      for (cT in 1:tabs) {
        for (cTT in 1:4) {
          firstlinespacing <- sprintf("%s ", firstlinespacing)
        }
      }
    }
    if ((!is.null(spaces)) & (spaces != 0)) {
      for (cT in 1:spaces) {
        firstlinespacing <- sprintf("%s ", firstlinespacing)
      }
    }
    if ((!is.null(indent)) & (indent > 0)) {
      # Add the indent to the text string
      for (cT in 1:indent) {
        firstlinespacing <- sprintf("%s ", firstlinespacing)
      }
    }

    # Create subsequent line spacing
    nextlinespacing <- ""
    if ((!is.null(tabs)) & (tabs != 0)) {
      for (cT in 1:tabs) {
        for (cTT in 1:4) {
          nextlinespacing <- sprintf("%s ", nextlinespacing)
        }
      }
    }
    if ((!is.null(spaces)) & (spaces != 0)) {
      for (cT in 1:spaces) {
        nextlinespacing <- sprintf("%s ", nextlinespacing)
      }
    }
    if ((!is.null(indent)) & (indent < 0)) {
      # Add the indent to the text string
      for (cT in 1:abs(indent)) {
        nextlinespacing <- sprintf("%s ", nextlinespacing)
      }
    }

    outputmat <- matrix(NA, nrow = 1, ncol = 1)
    outputvect <- matrix(NA, nrow = 1, ncol = 1)
    # Split the text based upon spaces
    tempsubstrings <- unlist(strsplit(as.character(textstring), " "))

    if (length(tempsubstrings) > 1) {
      rowindex <- 1
      startindex <- 1
      substringIndex <- 2
      while (substringIndex <= length(tempsubstrings)) {
        outputvect <- matrix(NA, nrow = 1, ncol = 1)
        spanstring <- paste(unlist(tempsubstrings[startindex:substringIndex]),collapse=" ")
        if (rowindex == 1) {
          spanstring <- sprintf("%s%s", firstlinespacing, spanstring)
        } else {
          spanstring <- sprintf("%s%s", nextlinespacing, spanstring)
        }
        if ((nchar(spanstring) > characters) | (nchar(spanstring) > wid) | (substringIndex==length(tempsubstrings))) {
         # if the character span exceeds the requested length or
         # if the character span exceeds the width of the screen or
         # if the end of the row has been reached
          outputvect[1,1] <- sprintf("%s",spanstring)
          outputmat <- base:cbind(outputmat, outputvect)
          rowindex <- rowindex + 1
          startindex <- substringIndex + 1
        }
        substringIndex <- substringIndex + 1
      }
    } else {
      spanstring <- sprintf("%s%s", firstlinespacing, tempsubstrings)
      if ((nchar(spanstring) < characters) & (nchar(spanstring) < wid)) {
        # text does not have any spaces and is shorter than the requested length
        outputvect[1,1] <- sprintf("%s",spanstring)
        outputmat <- base:cbind(outputmat, outputvect)
      } else {
        # text is too long so need to force it to wrap
        rowindex <- 1
        startindex <- 1
        substringIndex <- 2
        while (substringIndex <= nchar(tempsubstrings)) {
          outputvect <- matrix(NA, nrow = 1, ncol = 1)
          spanstring <- substr(tempsubstrings, startindex, substringIndex)
          if (rowindex == 1) {
            spanstring <- sprintf("%s%s", firstlinespacing, spanstring)
          } else {
            spanstring <- sprintf("%s%s", nextlinespacing, spanstring)
          }
          if ((nchar(spanstring) > characters) | (nchar(spanstring) > wid) | (substringIndex==nchar(tempsubstrings))) {
            # if the character span exceeds the requested length or
            # if the character span exceeds the width of the screen or
            # if the end of the row has been reached
            outputvect[1,1] <- sprintf("%s",spanstring)
            outputmat <- base:cbind(outputmat, outputvect)
            rowindex <- rowindex + 1
            startindex <- substringIndex + 1
          }
          substringIndex <- substringIndex + 1
        }
      }
    }

    for (cT in 2:ncol(outputmat)) {
      if (cT > 2) {
        cat(sprintf("\n"))
      }
      cat(outputmat[1,cT])
    }
    cat(sprintf("\n"))
  }
}
