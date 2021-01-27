
DGP_rlasso <- function(n, p, px, intercept = 1){
  
  X <- matrix(rnorm(n*p), ncol=p)
  colnames(X) <- paste("x", 1:p, sep="")
  beta <- c(rep(2,px), rep(0,p-px))
  y <- intercept + X %*% beta + rnorm(n)
  
  list(X = X, y = y, beta = beta)
}

DPG_lassoShooting <- function(n, p, px, lambda0 = 110, min = 0.85, max = 1.15){
  
  ret <- DGP_rlasso(n, p, px)
  
  loadings <- runif(p, min = min, max = max)
  lambda <- lambda0 * loadings
  
  list(X = ret$X, y = ret$y, beta = ret$beta, lambda = lambda, lambda0 = lambda0, loadings = loadings)
}

