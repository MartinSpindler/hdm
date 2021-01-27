context('Test rlasso vs manual')
library(hdm)
library(testthat)

set.seed(1)
sim_data <- DGP_rlasso_hdm_paper_sec4(n=5000, p=20)

test_cases <- expand.grid(post = c(TRUE, FALSE),
                          method = c('partialling out'))
test_cases['test_name'] <- apply(test_cases, 1, paste, collapse="_")

manual_partial_out <- function(data, x_cols, post){
  fmla.y <- as.formula(paste('y ~ ',paste(x_cols, collapse='+')))
  fmla.d <- as.formula(paste('d ~ ',paste(x_cols, collapse='+')))
  rY <- rlasso(fmla.y, data=data, post=post)$res
  rD <- rlasso(fmla.d, data=data, post=post)$res
  partial.fit <- lm(rY~rD)
  return(partial.fit)
}

patrick::with_parameters_test_that("Unit tests for rlassoEffect vs manual:", .cases = test_cases, {
  xnames <- colnames(sim_data)[-c(1,2)]
  
  partial.fit <- manual_partial_out(sim_data, xnames, post=post)
  
  Eff = rlassoEffect(as.matrix(sim_data[, xnames]), sim_data['y'], sim_data['d'],
                     method=method, post=post)
  
  expect_equal(summary(partial.fit)$coef['rD', 1:2], summary(Eff)$coef[,1:2], tolerance = 1e-8)
  
}
)
