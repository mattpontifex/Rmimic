context("descriptives")

test_that("descriptives works", {
  
  desc <- descriptives(data=mtcars, groupvariable=c("cyl"), verbose=FALSE)
  
  expect_true(as.integer(desc$Median[1]) == 26)
  expect_true(desc$Min[2] == 17.8)
  expect_true(is.na(desc$Frequencies[3]))
  
})