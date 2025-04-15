test_that("vanilla runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 family = fam)))
})

test_that("Component wise selection runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 varsel = TRUE, family = fam)))
})

test_that("Hierarchical selection runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 varsel = TRUE, groups = c(1,1,2,2,2), family = fam)))
})

test_that("Knots runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  suppressWarnings(knots10 <- fields::cover.design(dat$Z, nd = 10)$design)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 knots = knots10, family = fam)))
})

test_that("Rand intercept runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 id = rep(1:50, each = 2), family = fam)))
})

test_that("est.h runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 est.h = TRUE, family = fam)))
})

test_that("multiple groups runs", {
  set.seed(6)
  fam <- binomial(link = 'probit')
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                 X = dat$X, iter = 50, verbose = FALSE,
                                 K = 2, family = fam)))
})
