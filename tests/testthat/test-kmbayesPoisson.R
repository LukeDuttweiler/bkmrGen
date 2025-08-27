test_that("vanilla runs", {
  set.seed(6)
  fam <- poisson()
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                                  X = dat$X, iter = 50, verbose = FALSE,
                                                  family = fam)))
})

test_that("Component wise selection runs", {
  set.seed(6)
  fam <- poisson()
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                                  X = dat$X, iter = 50, verbose = FALSE,
                                                  varsel = TRUE, family = fam)))
})

test_that("multiple groups runs", {
  set.seed(6)
  fam <- poisson()
  dat <- SimData(family = fam)
  expect_no_error(suppressWarnings(tst <- kmbayes(y = dat$y, Z = dat$Z,
                                                  X = dat$X, iter = 50, verbose = FALSE,
                                                  K = 2, family = fam)))
})
