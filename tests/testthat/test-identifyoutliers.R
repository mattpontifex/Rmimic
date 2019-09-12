context("identifyoutliers")

test_that("identifyoutliers works", {
  
  testvect <- stats::runif(20,1,10)
  testvect[6] <- -1000
  testvect[10] <- 1000
  testvect <- identifyoutliers(testvect, iqrlimit = 3, verbose=FALSE)
  expect_true(is.na(testvect[6]))
  
})