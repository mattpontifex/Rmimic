context("RmimicAnova")

test_that("RmimicAnova works", {

  testout <- RmimicAnova(
                   data = PlantGrowth,
                   dependentvariable = "weight",
                   subjectid = NULL,
                   between = "group",
                   within = NULL, 
                   verbose=FALSE)
  
  expect_true(names(testout) == 4)
  expect_true(round(testout$stats$p[1], digits=4) == 0.0159)
  expect_true(round(testout$posthocttest$p[1], digits=4) == 0.249)
  
})