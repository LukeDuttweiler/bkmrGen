test_that("ExtractEsts works on gaussian", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2)

  expect_no_error(ests <- ExtractEsts(fitkm))
})


test_that("ExtractPIPs works on gaussian", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2)

  expect_no_error(ExtractPIPs(fitkm))
})
