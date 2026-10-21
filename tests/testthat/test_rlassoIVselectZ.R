context("Test rlassoIVselectZ")
library(hdm)
library(testthat)

# IV model without controls; y and d have nonzero means, true effect is -2
set.seed(1)
n <- 2000
z <- cbind(z1 = rnorm(n), z2 = rnorm(n))
u <- rnorm(n)
d <- 5 + 0.8*z[, 1] + 0.5*z[, 2] + 0.5*u + rnorm(n)
y <- 10 - 2*d + u


test_that("rlassoIVselectZ - second stage includes an intercept", {
  fit <- rlassoIV(x = NULL, d = d, y = y, z = z, select.X = FALSE, select.Z = TRUE)
  # both instruments are selected, so the estimator is 2SLS with an intercept
  expect_true(all(fit$selection.matrix))
  ref <- tsls(x = NULL, d = d, y = y, z = z, intercept = TRUE, homoscedastic = FALSE)
  expect_equal(unname(fit$coefficients), unname(ref$coefficients[1, 1]), tolerance = 1e-8)
  expect_equal(unname(fit$se), unname(ref$se[1]), tolerance = 1e-8)
  expect_equal(unname(fit$coefficients), -2, tolerance = 0.1)
})

test_that("rlassoIVselectZ - estimate does not depend on the location of y", {
  fit <- rlassoIVselectZ(x = NULL, d = d, y = y, z = z)
  fit.shift <- rlassoIVselectZ(x = NULL, d = d, y = y + 100, z = z)
  expect_equal(fit.shift$coefficients, fit$coefficients, tolerance = 1e-8)
  expect_equal(fit.shift$se, fit$se, tolerance = 1e-8)
})

test_that("rlassoIVselectZ - robust variance when controls are left out of the first stage", {
  set.seed(3)
  n <- 500
  p <- 10
  zc <- cbind(z1 = rnorm(n), z2 = rnorm(n))
  xc <- matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  uc <- rnorm(n)
  dc <- zc[, 1] + 0.5*zc[, 2] + xc[, 1] + 0.15*rowSums(xc[, -1]) + uc + rnorm(n)
  yc <- -2*dc + rowSums(xc) + uc
  fit <- rlassoIV(x = xc, d = dc, y = yc, z = zc, select.X = FALSE, select.Z = TRUE)
  # some controls are not selected, so t(Dhat) %*% d is not symmetric
  expect_true(any(!fit$selection.matrix[colnames(xc), 1]))
  # the estimator is IV of y on (d, x) with instruments (first-stage fit, x)
  dhat <- predict(rlasso(x = cbind(zc, xc), y = dc, post = TRUE, intercept = TRUE))
  ref <- tsls(x = xc, d = dc, y = yc, z = dhat, intercept = TRUE, homoscedastic = FALSE)
  expect_equal(unname(fit$coefficients), unname(ref$coefficients[1, 1]), tolerance = 1e-8)
  expect_equal(unname(fit$se), unname(ref$se[1]), tolerance = 1e-8)
})

test_that("rlassoIVselectZ - exact solve with ill-conditioned controls", {
  data(EminentDomain)
  ze <- EminentDomain$logGDP$z
  xe <- EminentDomain$logGDP$x
  ye <- EminentDomain$logGDP$y
  de <- EminentDomain$logGDP$d
  fit <- rlassoIV(x = xe, d = de, y = ye, z = ze, select.X = FALSE, select.Z = TRUE)
  # Frisch-Waugh: IV of y on d with instrument dhat, after partialling out x
  # (x contains a constant column)
  dhat <- predict(rlasso(x = cbind(ze, xe), y = de, post = TRUE, intercept = TRUE))
  qx <- qr(xe)
  expect_equal(qx$rank, ncol(xe))
  dh.r <- qr.resid(qx, dhat)
  d.r <- qr.resid(qx, de)
  y.r <- qr.resid(qx, ye)
  b <- sum(dh.r*y.r)/sum(dh.r*d.r)
  e <- y.r - b*d.r
  expect_equal(unname(fit$coefficients), b, tolerance = 1e-6)
  expect_equal(unname(fit$se), sqrt(sum(dh.r^2*e^2))/abs(sum(dh.r*d.r)), tolerance = 1e-6)
  # the estimate does not depend on the units of the controls
  sc <- 10^seq(-3, 3, length.out = ncol(xe))
  fit.sc <- rlassoIV(x = sweep(xe, 2, sc, "*"), d = de, y = ye, z = ze,
                     select.X = FALSE, select.Z = TRUE)
  expect_equal(fit.sc$coefficients, fit$coefficients, tolerance = 1e-6)
  expect_equal(fit.sc$se, fit$se, tolerance = 1e-6)
})
