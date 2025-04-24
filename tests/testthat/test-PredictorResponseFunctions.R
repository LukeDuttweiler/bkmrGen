test_that("Predictor response functions work for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
  expect_no_error(pred.resp.univar <- PredictorResponseUnivar(fit = fitkm))
  expect_no_error(pred.resp.bivar12 <- PredictorResponseBivarPair(fit = fitkm, min.plot.dist = 1, ngrid = 10))
  expect_no_error(pred.resp.bivar <- PredictorResponseBivar(fit = fitkm, min.plot.dist = 1, ngrid = 10))
  expect_no_error(pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.df = pred.resp.bivar,
                                                                         Z = Z, qs = c(0.1, 0.5, 0.9)))
})
