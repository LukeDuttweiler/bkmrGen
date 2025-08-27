test_that("ExtractEsts works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2)

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with est.h = TRUE", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                   est.h = TRUE)

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with probit", {
  set.seed(111)
  fam <- binomial(link = 'probit')
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                                    family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with logit", {
  set.seed(111)
  fam <- binomial(link = 'logit')
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                                    family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with poisson", {
  set.seed(111)
  fam <- poisson()
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                                    family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with est.h = TRUE, probit", {
  set.seed(111)
  fam <- binomial(link = 'probit')
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                                    est.h = TRUE, family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with est.h = TRUE, logit", {
  set.seed(111)
  fam <- binomial(link = 'logit')
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                                    est.h = TRUE, family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})

test_that("ExtractEsts works with est.h = TRUE, poisson", {
  set.seed(111)
  fam <- poisson()
  dat <- SimData(n = 50, M = 4, family = fam)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- suppressWarnings(kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                   est.h = TRUE, family = fam))

  expect_no_error(ests <- ExtractEsts(fitkm))
})


test_that("ExtractPIPs works for default", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2)

  expect_no_error(ExtractPIPs(fitkm))
})

test_that("ExtractPIPs works with est.h = TRUE", {
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE, K = 2,
                   est.h = TRUE)

  expect_no_error(ExtractPIPs(fitkm))
})
