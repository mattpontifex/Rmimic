tryCatch(library(testthat), error=function(e){install.packages("testthat"); library(testthat)}) # test_check
test_check("Rmimic")