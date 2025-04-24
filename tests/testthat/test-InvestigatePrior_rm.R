test_that("InvestigatePrior and Plot works", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  expect_no_error(priorfits <- InvestigatePrior(y = y, Z = Z, X = X, q.seq = c(2, 1/2, 1/4, 1/16)))
  expect_no_error(PlotPriorFits(y = y, Z = Z, X = X, fits = priorfits))
})
