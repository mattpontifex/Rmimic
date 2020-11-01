context("RmimicMLAnova")

test_that("RmimicMLAnova works", {
  
  testout <- RmimicMLAnova(data = alertness[which(alertness$Condition == 'Condition2'),],
       dependentvariable = 'Alertness',
       subjectid = 'PartID', 
       between = 'Group', 
       within = c('Time'), 
       randomintercept = c('PartID'),
       posthoc = 'False Discovery Rate Control', 
       planned = NULL, 
       studywiseAlpha = 0.05, confidenceinterval = 0.95, verbose=FALSE)
  
  
  expect_true(round(testout$stats$p.value[1], digits=4) == 0.0585)
  expect_true(round(testout$posthocttest$p.value[1], digits=4) == 0.0055)
  
})