context("RmimicLMcontrast")

test_that("RmimicLMcontrast works", {

  basefit <- lm(mpg ~ am + wt, data = mtcars)
  fit <- lm(mpg ~ am + wt + qsec, data = mtcars)
  regresult <- RmimicLMcontrast(basefit, fit, confidenceinterval=0.95, studywiseAlpha=0.05, verbose=FALSE)
  
  expect_true(length(regresult) == 3)
  expect_true(regresult$stats$DFd[1] == 29)
  
})