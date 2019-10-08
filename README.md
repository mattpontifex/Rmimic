# Rmimic
An R package with miscellaneous R functions that are useful to mimic functionalities of popular statistics software packages. This package deviates from the typical R ethos such that most functions will output results in the console window (results are also available as environment variables; and output can be suppressed using verbose calls to the function).

To use this package, from R run the following commands:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*devtools::install_github("mattpontifex/Rmimic"); library(Rmimic)*

# Function List
**ttest2text**: Function that takes a t-test result from the stats package and outputs the t-test result for use in an APA style manuscript (i.e., t(21) = 2.1, p = 0.001) with proper rounding.
**descriptives**: Function that computes SPSS style descriptive statistics and frequencies.
**clusterthreshold1d**: Function that calculates contiguous clusters of locations in a 1D array that are above or below some threshold and of some minimum cluster size.
**multipleimputation**: Function that uses the mice package to replace missing data points.
