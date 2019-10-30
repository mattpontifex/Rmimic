context("RmimicAnova")

test_that("RmimicAnova works", {

  testout <- RmimicAnova(
                   data = PlantGrowth,
                   dependentvariable = "weight",
                   subjectid = NULL,
                   between = "group",
                   within = NULL, 
                   verbose=FALSE)
  
  #expect_true(length(testout) > 1)
  #expect_true(testout$descriptives$N[1] == 10)
  #expect_true(testout$stats$df[1] == 18)
  
})