test_that("SamplePred works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)

  med_vals <- apply(Z, 2, median)
  Znew <- matrix(med_vals, nrow = 1)
  h_true <- dat$HFun(Znew)
  set.seed(111)
  expect_no_error(samps3 <- SamplePred(fitkm, Znew = Znew, Xnew = cbind(0)))
})
