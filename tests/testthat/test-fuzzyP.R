context("fuzzyP")

test_that("fuzzyP works", {
  
  outPvalue <- fuzzyP(0.04)
  
  expect_equal(outPvalue$interpret, 0.04)
  expect_equal(outPvalue$report, "0.04")
  expect_equal(outPvalue$modifier, "=")
  
})