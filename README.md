# Rmimic
An R package with miscellaneous R functions that are useful to me.

To use this package, from R run the following commands:
tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})
tryCatch(library(Rmimic), error=function(e){devtools::install_github("mattpontifex/Rmimic"); library(Rmimic)})
