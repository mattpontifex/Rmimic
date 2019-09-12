context("determineallpossiblecombinations")

test_that("determineallpossiblecombinations works", {
  
  outlist <- determineallpossiblecombinations(c('A', 'B'))
  
  expect_equal(outlist[1], "A")
  expect_equal(outlist[2], "B")
  expect_equal(outlist[3], "A:B")
  
})