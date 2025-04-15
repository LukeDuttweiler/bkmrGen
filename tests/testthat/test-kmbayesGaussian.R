test_that("vanilla runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z, X = dat$X, iter = 50, verbose = FALSE))
})

test_that("Component wise selection runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 varsel = TRUE))
})

test_that("Hierarchical selection runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 varsel = TRUE, groups = c(1,1,2,2,2)))
})

test_that("Knots runs", {
  dat <- SimData()
  suppressWarnings(knots10 <- fields::cover.design(dat$Z, nd = 10)$design)
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 knots = knots10))
})

test_that("Rand intercept runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 id = rep(1:50, each = 2)))
})

test_that("est.h runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 est.h = TRUE))
})

test_that("multiple groups runs", {
  dat <- SimData()
  expect_no_error(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 K = 2))
})
