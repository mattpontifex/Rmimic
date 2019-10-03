# Rmimic
An R package with miscellaneous R functions that are useful to mimic functionalities of popular statistics software packages. This package deviates from the typical R mythos such that most functions will output results in the console window (results are also available as environment variables; and output can be suppressed using verbose calls to the function).

To use this package, from R run the following commands:

tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})

devtools::install_github("mattpontifex/Rmimic"); library(Rmimic)
