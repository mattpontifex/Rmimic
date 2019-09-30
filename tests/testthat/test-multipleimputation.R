context("multipleimputation")

test_that("multipleimputation works", {
  
  data <- multipleimputation(airquality, imputations=10)
  
  expect_true(!is.na(data$Ozone[5]))
  expect_true(!is.na(data$Ozone[10]))
  
})