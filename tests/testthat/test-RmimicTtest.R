context("RmimicTtest")

test_that("RmimicTtest works", {
  
  testout <- RmimicTtest(PlantGrowth, 
                         dependentvariable='weight', 
                         subjectid=NULL, 
                         between='group', 
                         within=NULL, 
                         nonparametric=FALSE,
                         posthoc="Holm-Bonferroni", 
                         verbose=FALSE)
  
  expect_true(length(testout) > 1)
  expect_true(testout$descriptives$N[1] == 10)
  expect_true(testout$stats$df[1] == 18)
  
})