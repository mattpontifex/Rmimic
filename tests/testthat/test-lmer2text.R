context("lmer2text")

test_that("lmer2text works", {
  
  mockdatabase <- PlantGrowth
  mockdatabase$observation <- rep_len(c(1,2),30)
  mockdatabase$PartID <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15)
  
  pkgcond::suppress_conditions(suppressWarnings(fit <- lmerTest::lmer(weight ~ group*observation + (1 | PartID), data = mockdatabase)))
  result <- lmer2text(fit, df="Kenward-Roger", numparticipants=15, numfactors=3)
  
  expect_true(length(result) > 1)
  expect_true(result$ANOVA$DFn[1] == 2)
  
})