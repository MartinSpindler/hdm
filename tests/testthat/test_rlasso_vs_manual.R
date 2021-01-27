context('Test rlasso vs manual')
library(hdm)
library(testthat)

set.seed(1)
sim_data <- list()
sim_data[[1]] <- DGP_rlasso_hdm_paper_sec4(n=5000, p=20)

xx <- DGP_rlasso(200, 100, 10)
colnames(xx$X) <- c("d",paste("x",1:(100-1),sep=""))
sim_data[[2]] <- data.frame(y= xx$y, xx$X)

test_cases <- expand.grid(post = c(TRUE, FALSE),
                          method = c('partialling out', 'double selection'),
                          i_data = 1:(length(sim_data)))

test_cases['test_name'] <- apply(test_cases, 1, paste, collapse="_")

manual_partial_out <- function(data, x_cols, post){
  fmla.y <- as.formula(paste('y ~ ',paste(x_cols, collapse='+')))
  fmla.d <- as.formula(paste('d ~ ',paste(x_cols, collapse='+')))
  rY <- rlasso(fmla.y, data=data, post=post)$res
  rD <- rlasso(fmla.d, data=data, post=post)$res
  partial.fit <- lm(rY~rD)
  return(partial.fit)
}

manual_double_selection <- function(data, x_cols, post){
  fmla.y <- as.formula(paste('y ~ ',paste(x_cols, collapse='+')))
  fmla.d <- as.formula(paste('d ~ ',paste(x_cols, collapse='+')))
  betaY <- rlasso(fmla.y, data=data, post=post)$coefficients
  betaD <- rlasso(fmla.d, data=data, post=post)$coefficients
  
  sel_x = x_cols[abs(betaY[x_cols])>0 | abs(betaD[x_cols])>0]
  fmla.double.sel <- as.formula(paste('y ~ ', paste(c("d", sel_x), collapse='+')))
  partial.fit <- lm(fmla.double.sel, data)
  return(partial.fit)
}

patrick::with_parameters_test_that("Unit tests for rlassoEffect vs manual:", .cases = test_cases, {
  this_data <- sim_data[[i_data]]
  
  xnames <- colnames(this_data)[-c(1,2)]
  
  if (method == 'partialling out') {
    partial.fit <- manual_partial_out(this_data, xnames, post=post)
  } else if (method == 'double selection') {
    partial.fit <- manual_double_selection(this_data, xnames, post=post)
  }
  
  Eff = rlassoEffect(as.matrix(this_data[, xnames]), this_data['y'], this_data['d'],
                     method=method, post=post)
  
  if (method == 'partialling out') {
    expect_equal(summary(partial.fit)$coef['rD', 1:2],
                 summary(Eff)$coef[,1:2], tolerance = 1e-8)
  } else if (method == 'double selection') {
    expect_equal(partial.fit$coefficients['d'],
                 Eff$alpha, tolerance = 1e-8)
  }
}
)
