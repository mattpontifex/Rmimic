context("clusterthreshold1d")

test_that("clusterthreshold1d works", {
  
  testvect <- stats::runif(20,1,10)
  testvect[4:10] <- c(100, 100, 100, 100, 100, 100, 100)
  testvect <- clusterthreshold1d(testvect, crit = 50, clustersize = 5, direction = 'GreaterThan')
  
  expect_equal(testvect[1], 0)
  expect_equal(testvect[20], 0)
  expect_equal(testvect[6], 1)
  expect_equal(testvect[8], 1)
  
})