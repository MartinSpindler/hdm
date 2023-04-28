context("AutoDML")
library(hdm)
library(testthat)
library(data.table)

theta = 1

DGP_bintreat = function(n_obs = 500, dim_x = 20, theta = 0, R2_d = 0.5,
  R2_y = 0.5) {
  v = runif(n_obs)
  zeta = rnorm(n_obs)
  cov_mat = toeplitz(0.5^(0:(dim_x - 1)))
  x = matrix(rnorm(n_obs * dim_x), nrow = n_obs, ncol = dim_x)
  beta = 1 / (1:dim_x)^2
  b_sigma_b = beta %*% cov_mat %*% beta
  c_y = c(R2_y / ((1 - R2_y) * b_sigma_b))
  c_d = c(pi^2 / 3 * R2_d / ((1 - R2_d) * b_sigma_b))
  xx = exp(x %*% (beta * c_d))
  d = 1 * ((xx / (1 + xx)) > v)
  y = d * theta + d * x %*% (beta * c_y) + zeta
  colnames(x) = paste0("X", 1:dim_x)
  colnames(y) = "Y"
  colnames(d) = "D"
  data = data.table(x, y, d)
  return(data)
}


set.seed(2)
df <- DGP_bintreat(100, 10, theta)

# Result pre-refactoring
dictionary <- function(d, z) {
  return(c(1, d, z))
}

Xnames = colnames(df)[grep("X", colnames(df))]
backend_binD = DataATEAutoDML(x = Xnames, d = "D", y = "Y",
                              data = df, dict = dictionary)

set.seed(2)
auto_ate = rlassoAutoDML(backend_binD, prelim_est = FALSE, est_type = "ATE")
set.seed(2)
auto_ate_prelim = rlassoAutoDML(backend_binD, prelim_est = TRUE, est_type = "ATE")


# df = data.frame(df)
# Xs = as.matrix(df[grep("X", names(df))])
# D = df[, "D"]
# Y = df[, "Y"]
# ate_autodml = rlassoATEAutoDML(Xs, D, Y, dictionary)
# > summary(ate_autodml)
# Estimation and significance testing of the treatment effect
# Type: ATE
#    coeff.    se. t-value  p-value
# TE 1.3077 0.2484   5.264 1.41e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#####

# TODO: Fix tests below!


test_that("rlassoEffect - Input check x, y and d", {
  expect_is(rlassoEffect(X, y, d = X[, 1]), "rlassoEffects")
  expect_is(rlassoEffect(X, as.vector(y), d = X[, 1]), "rlassoEffects")
  expect_is(rlassoEffect(X[, 1, drop = FALSE], y, d = X[, 1]), "rlassoEffects")
  expect_is(rlassoEffect(X[, 1, drop = FALSE], as.vector(y), d = X[, 1]), "rlassoEffects")
})

test_that("rlassoEffect - Input check I3", {
  expect_is(rlassoEffect(X, y, d = X[, 1], I3 = c(rep(TRUE, 2), rep(FALSE, 2), TRUE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], I3 = c(rep(TRUE, 55), rep(FALSE, 44), TRUE)), "rlassoEffects")
})

test_that("rlassoEffect - Input check post, intercept and normalize", {
  expect_is(rlassoEffect(X, y, d = X[, 1], post = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], intercept = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], normalize = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], normalize = FALSE, intercept = FALSE), "rlassoEffects")
})

test_that("rlassoEffect - Input check penalty", {
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = TRUE, X.dependent.lambda = FALSE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = FALSE, X.dependent.lambda = TRUE, numSim = 4000)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = TRUE, X.dependent.lambda = TRUE, numSim = 4000)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = "none", X.dependent.lambda = FALSE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(homoscedastic = "none", X.dependent.lambda = TRUE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], intercept = FALSE, penalty = list(homoscedastic = "none", X.dependent.lambda = FALSE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], penalty = list(
    homoscedastic = FALSE, X.dependent.lambda = FALSE,
    lambda.start = NULL, c = 1.1, gamma = 0.1
  )), "rlassoEffects")
})

test_that("rlassoEffect - Input check control", {
  expect_is(rlassoEffect(X, y, d = X[, 1], control = list(numIter = 15, tol = 10^-4, threshold = 10^-3)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], control = list(numIter = 25)), "rlassoEffects")
})
