context("Test LassoShooting.fit")
library(hdm)
library(testthat)

set.seed(2)
ret <- DPG_lassoShooting(200, 100, 10, 110)
X <- ret$X
y <- ret$y
beta <- ret$beta
lambda <- ret$lambda
rm(ret)


test_that("LassoShooting.fit - Input check x, y and lambda",{
  expect_is(LassoShooting.fit(X, y, lambda), "list")
  expect_is(LassoShooting.fit(X, as.vector(y), lambda), "list")
  expect_is(LassoShooting.fit(X[, 1, drop = FALSE], y, lambda), "list")
  expect_is(LassoShooting.fit(X[, 1, drop = FALSE], as.vector(y), lambda), "list")
})


test_that("LassoShooting.fit - Input check control, XX, Xy and beta start",{
  expect_is(LassoShooting.fit(X, y, lambda, control = list(maxIter = 150, optTol = 10^(-4), zeroThreshold = 10^(-5))), "list")
  expect_is(LassoShooting.fit(X, y, lambda, XX = (t(X) %*% X) * 0.8), "list")
  expect_is(LassoShooting.fit(X, y, lambda, Xy = (t(X) %*% y) * 0.8), "list")
  expect_is(LassoShooting.fit(X, y, lambda, beta.start = rep(1,100)), "list")
})