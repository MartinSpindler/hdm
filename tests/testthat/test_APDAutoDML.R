context("AutoDML")
library(hdm)
library(testthat)
library(data.table)
library(MASS)

DGP_contD <- function(n, p){
  X <- matrix(rnorm(n * p), ncol = p)
  DT <- data.table(X)
  theta = 1/(1:p)^2 
  # creating the treatment 
  DT$D <- .1*(X%*%theta) + rbeta(n, 1, 7)
  DT$y <- DT$D  + DT$D^2 + DT$D^3 + DT$D*DT$V1 +  .1*(X%*%theta) + rnorm(n, 0,1)
  return(list("data" = DT, "true" = mean(1 + 2*DT$D + 3*DT$D^2 + DT$V1)))
}

DGP_panel <- function(N, N_T, p){
    a =rnorm(N, 1, 1)
    fixed =rep(1:N, each = N_T)
    time =rep(1:N_T,  N)
    DT <- data.table(fixed, time, a = rep(a, each = N_T))
    size=p
    k = diag(size) 
    sigma = k
    theta = 1/(1:p)^2
    X = (Reduce(rbind, lapply(1:nrow(DT), function(i){mvrnorm(n =1, rep(DT$a[i], size), sigma)})))
    DT <- cbind(DT, setNames(data.table(X), paste0("X_", 1:p)))
    DT$D <- .1*(X%*%theta) + rbeta(N, 1, 7)
    DT$y <- DT$a + DT$D  + DT$D^2 + DT$D^3 + DT$D*DT$X_1 +  .1*(X%*%theta) + rnorm(N, 0,1)
    
    return(list("data" = DT, "true" = mean(1 + 2*DT$D + 3*DT$D^2 + DT$X_1)))
}

test_cases = expand.grid(
  panel_data = c(FALSE, TRUE),
  L = c(2, 5),
  gamma_learner = c("rlasso", "cv.glmnet"),
  est_type = "APD",
  poly_degree = c(2, 4),
  interactions = c(FALSE, TRUE),
  # weights = c(NULL, weights_vec),
  debiased = c(TRUE, FALSE),
  D_LB = c(0, 0.01),
  D_add = c(0.2, 0.25),
  stringsAsFactors = FALSE)

test_cases[".test_name"] = apply(test_cases, 1, paste, collapse = "_")

patrick::with_parameters_test_that("Unit tests for rlassoAutoDML for ATE:",
  .cases = test_cases, {
    
    set.seed(2)
    if (!panel_data) {
      data_gen <- DGP_contD(500, 10)
      df <- data_gen$data
      df$intercept <- rep(1, nrow(df))
      Xnames <- colnames(df)[grep("V", colnames(df))]
      backend_contD = DataAPDAutoDML(x = Xnames, d = "D", y = "y", data = df,
                                     poly = poly_degree, interactions = interactions)
    } else {
      dgp_panel <- DGP_panel(1000, 10, 4)
      DT_panel <- dgp_panel$data
      Xnames = colnames(DT_panel)[grep("X", colnames(DT_panel))]
      backend_contD <- DataAPDAutoDML(x = Xnames, d = "D", y = "y",
                                            data = DT_panel,
                                            poly_order = poly_degree,
                                            interactions = interactions,
                                            unit = "fixed",
                                            time = "time")
    }
    
    auto_apd = rlassoAutoDML(backend_contD, L = L,
                             gamma_learner = gamma_learner,
                             est_type = est_type,
                             # weights = weights,
                             debiased = debiased,
                             D_LB = D_LB,
                             D_add = D_add)
    
    
    expect_is(backend_contD, "DataAPDAutoDML")
    expect_is(auto_apd, "rlassoAutoDML")
    expect_is(auto_apd$te, "numeric")
  }
)
