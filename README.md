# Rmimic
An R package with miscellaneous R functions that are useful to mimic functionalities of popular statistics software packages. This package deviates from the typical R ethos such that most functions will output results in the console window (results are also available as environment variables; and output can be suppressed using verbose calls to the function).

To use this package, from R run the following commands:

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;devtools::install_github("mattpontifex/Rmimic"); library(Rmimic)

# Main Function List
* **clusterthreshold1d**: Function that calculates contiguous clusters of locations in a 1D array that are above or below some threshold and of some minimum cluster size (i.e., a cluster of 30 points all below 0.05).
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;tempdata <- c(0.2, 0.3, 0.05, 0.04, 0.06, 0.08, 0.009, 0.05, 0.02, 0.03, 0.08, 0.1, 0.4)  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;clusterthreshold1d(tempdata, crit = 0.05, clustersize = 3, direction = 'LessThan')  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;returns: *c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)*

* **descriptives**: Function that computes SPSS style descriptive statistics and frequencies.
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;tempdata <- data.frame("Group" = sample(1:2,100, replace=TRUE), "X" = runif(100), "Y" = runif(100))  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;desc <- descriptives(tempdata, groupvariable=c("Group"), verbose=TRUE)  

* **identifyoutliers**: Function to identify outliers based upon the interquartile range (as SPSS does for boxplots) and replace those values with NA.
* **multipleimputation**: Function that uses the mice package to replace missing data points.
* **ttest2text**: Function that takes a t-test result from the stats package and outputs the t-test result for use in an APA style manuscript (i.e., t(18) = 2.3, p = 0.031) with proper rounding.
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;comparison1 <- c(2, 3, 4, 6, 8, 9, 10, 11, 13, 15); comparison2 <- c( 5, 7, 9, 10, 13, 15, 16, 17, 18, 20)  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ttest <- stats::t.test(x= comparison1, y=comparison2, paired=FALSE, var.equal=TRUE)  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ttest$effectsize <- ttest$statistic * sqrt((1/length(comparison1)) + (1/length(comparison2)))  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;tempout <- ttest2text(ttest, verbose=TRUE)  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;returns: *t(18) = 2.3, p = 0.031, ds = -1.04.*


# Accessory Function List
* **decimalplaces**: Function to obtain the number of decimal places of precision in a vector.
* **determineallpossiblecombinations**: Function to determine all possible combinations of an input array. For instance, an array containing A, B, and C could be assessed looking at A, B, C, A:B, A:C, B:C, or A:B:C.
* **fuzzyP**: Function to round P values for reporting. Because there is no reason to report p = 0.912 to three digits of precision.
* **typewriter**: Function to control the text outputted to the console to create a consistent indent.
