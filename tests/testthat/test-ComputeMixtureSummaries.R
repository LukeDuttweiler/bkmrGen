test_that("OverallRiskSummaries works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)

  expect_no_error(risks.overall <- OverallRiskSummaries(fit = fitkm, qs = seq(0.25, 0.75, by = 0.05),
  q.fixed = 0.5, method = "exact"))
})

test_that("SingleVarRiskSummaries works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)

  expect_no_error(risks.singvar <- SingVarRiskSummaries(fit = fitkm, method = "exact"))
})

test_that("SingleVarIntSummaries works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)

  expect_no_error(SingVarIntSummaries(fit = fitkm, method = "exact"))
})
