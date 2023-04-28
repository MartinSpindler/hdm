context("AutoDML")
library(hdm)
library(testthat)
library(data.table)

DGP_contD <- function(n, p){
  X <- matrix(rnorm(n * p), ncol = p)
  DT <- data.table(X)
  theta = 1/(1:p)^2 
  # creating the treatment 
  DT$D <- .1*(X%*%theta) + rbeta(n, 1, 7)
  DT$y <- DT$D  + DT$D^2 + DT$D^3 + DT$D*DT$V1 +  .1*(X%*%theta) + rnorm(n, 0,1)
  return(list("data" = DT, "true" = mean(1 + 2*DT$D + 3*DT$D^2 + DT$V1)))
}


set.seed(2)
data_gen <- DGP_contD(100, 10)
df <- data_gen$data

# add intercept
df$intercept <- rep(1, nrow(df))
theta_true <- data_gen$true

poly_degree = 3

Xnames = colnames(df)[grep("V", colnames(df))][-1]
Xmanual = colnames(df)[grep("V", colnames(df))][1]

backend_contD = DataAPDAutoDML(x = Xnames, d = "D", y = "y", x_manual = Xmanual, data = df, poly = poly_degree, intercept = "intercept")

auto_apd = rlassoAutoDML(backend_contD, est_type = "APD")




# cases
poly = c(2, 3)
interactions = c(TRUE, FALSE)

#####



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
