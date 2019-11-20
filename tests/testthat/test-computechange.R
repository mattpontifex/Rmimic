context("computechange")

test_that("computechange works", {
  
  mockdatabase <- data.frame("ID" = rep_len(1:20,60), "Time" = c(rep_len("Time1",20),rep_len("Time2",20),rep_len("Time3",20)), "X" = runif(60))
  mockdatabase <- computechange(mockdatabase, dependentvariable='X', subjectid='ID', within='Time')
  expect_true(colnames(mockdatabase)[length(colnames(mockdatabase))] == 'X.changeinTimerelativetoTime1')
  
})