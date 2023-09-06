context("AutoDML (ATE)")
library(hdm)
library(testthat)
library(data.table)

theta = 1
n_obs = 500

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

# Use different dictionaries
dictionary1 <- function(d, z) {
  return(c(1, d, z))
}

dictionary2 <- function(d, z) {
  return(c(1, d, z, d*z))
}

dictionary3 <- function(d, z) {
  return(c(1, d, z, d*z, z*z))
}

set.seed(2)
df <- DGP_bintreat(n_obs, 10, theta)
weights_vec = rep(c(0,1),n_obs/2)

test_cases = expand.grid(
  dictionary_name = c("dict1"),#, "dict2", "dict3"),
  L = c(2, 5),
  gamma_learner = c("rlasso"),#, "cv.glmnet"),
  est_type = "ATE",
  prelim_est = c(TRUE, FALSE),
  # weights = c(NULL, weights_vec),
  debiased = c(TRUE), #, FALSE),
  D_LB = c(0),
  D_add = c(0.2),
  stringsAsFactors = FALSE)

test_cases[".test_name"] = apply(test_cases, 1, paste, collapse = "_")

patrick::with_parameters_test_that("Unit tests for rlassoAutoDML for ATE:",
  .cases = test_cases, {
    
    if (dictionary_name == "dict1") {
      dictionary = dictionary1
    }
    
    if (dictionary_name == "dict2") {
      dictionary = dictionary2
    }
    
    if (dictionary_name == "dict3") {
      dictionary = dictionary3
    }
    
  Xnames = colnames(df)[grep("X", colnames(df))]
  backend_binD = DataATEAutoDML(x = Xnames, d = "D", y = "Y",
                              data = df, dict = dictionary)

  auto_ate = rlassoAutoDML(backend_binD, L = L,
                           gamma_learner = gamma_learner,
                           est_type = est_type,
                           prelim_est = prelim_est,
                          # weights = weights,
                           debiased = debiased,
                           D_LB = D_LB,
                           D_add = D_add)
  
  
  expect_is(backend_binD, "DataATEAutoDML")
  expect_is(auto_ate, "rlassoAutoDML")
  expect_is(auto_ate$te, "numeric")

})
