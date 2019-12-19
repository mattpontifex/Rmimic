context("lmer2text")

test_that("lmer2text works", {

  fit <- lmerTest::lmer(Alertness ~ Group*Drug*Dose + (1 | PartID), data=elashoff)
  result <- lmer2text(fit, df="Kenward-Roger", numparticipants=16, numfactors=4)
  
  expect_true(length(result) > 1)
  expect_true(result$ANOVA$DFn[1] == 1)
  expect_true(result$ANOVA$p.value[2] < 0.05)
  
})