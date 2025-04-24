test_that("ComputePostmeanHnew works for gaussian", {
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
  expect_no_error(h_est1 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "approx"))
  expect_no_error(h_est2 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "exact"))
})
