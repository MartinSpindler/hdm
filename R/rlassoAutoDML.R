#' Auto DML based on rlasso
#' 
#' Implements the Double ML algorithm introduced in for estimating causal effects
#' This method was first introdued in https://arxiv.org/abs/1608.00060 which required
#' manual calculation of the Riesz representer. In Auto DML (from https://arxiv.org/abs/1809.05224)
#' the Riesz representer is estimated automatically and doesn't need to be explicitly 
#' computed.
#' 
#' This implementation 
#'
#' @param Y A vector of outputs
#' @param D A vector of treatment values
#' @param X A matrix of covariates
#' @param dict A dictionary
#' function of (d,z) that maps to a vector
#' default is (1,d,z)
#' @param bias debiased vs. biased results
#' @param D_LB Lower bound on D (default 0)
#' @param D_add value added to D (default 0.2)
#' @param L number of folds data is split into (default 5)
#' @param max_iter maximum iterations of Lasso (default 10)
#' @return list with average treatment effect and standard error
#' @examples
#' # data = simulate_data(500)
#' #
#' # Y = data[[1]]
#' # D = data[[2]]
#' # X = data[[3]]
#'
#' # rlassoAutoDML(Y,D,X,dict = b2)
#' # rlassoAutoDML(Y, T, X)
#' @export
#' @rdname rlassoDML
rlassoAutoDML <- function(Y, D, X, dict = NULL, D_LB = 0, D_add = 0.2,
                          bias = FALSE, L = 5, max_iter = 10) {
  
  if (is.null(dict)) {
    dict = function(d,z){
      return(c(1,d,z))
    }
    
  }
  p <- length(dict(D[1], X[1, ]))
  
  # p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  # TODO: is this comment true or is p0 = dim(X0)/4 used or can we drop this?
  p0 <- ceiling(p / 4)
  if (p > 60) {
    p0 <- ceiling(p / 40)
    
  }
  n <- nrow(X)
  folds <- split(sample(n, n, replace = FALSE), as.factor(1:L))
  
  Psi_tilde <- numeric(0)
  
  for (l in 1:L) {
    Y.l <- Y[folds[[l]]]
    Y.nl <- Y[-folds[[l]]]
    
    T.l <- D[folds[[l]]]
    T.nl <- D[-folds[[l]]]
    
    X.l <- X[folds[[l]], ]
    X.nl <- X[-folds[[l]], ]
    
    
    n.l <- length(T.l)
    n.nl <- length(T.nl)
    
    # get stage 1 (on nl)
    rho_hat <- RMD_stable(Y.nl, T.nl, X.nl, p0, D_LB, D_add, max_iter, dict)
    
    alpha_hat <- function(d, z) {
      return(dict(d, z) %*% rho_hat)
    }
    
    n <- nrow(X.nl)
    p <- length(dict(T.nl[1], X.nl[1, ]))
    # Apply the dictionary b to W
    B <- matrix(0, n, p)
    for (i in 1:n) {
      B[i, ] <- dict(T.nl[i], X.nl[i, ])
    }
    
    gamma_coeff <- rlasso(B, Y.nl, intercept = F)$coefficients
    gamma_hat <- function(d, z) {
      return(dict(d, z) %*% gamma_coeff)
    }
    
    print(paste0("fold: ", l))
    
    # get stage 2 (on l)
    # psi_star
    Psi_tilde.l <- rep(0, n.l)
    for (i in 1:n.l) {
      if (bias) { # plug-in
        Psi_tilde.l[i] <- psi_tilde_bias(Y.l[i], T.l[i], X.l[i, ], m, alpha_hat, gamma_hat) # without subtracting theta_hat
      } else { # DML
        Psi_tilde.l[i] <- psi_tilde(Y.l[i], T.l[i], X.l[i, ], m, alpha_hat, gamma_hat) # without subtracting theta_hat
      }
    }
    Psi_tilde <- c(Psi_tilde, Psi_tilde.l)
  }
  
  # point estimation
  ate <- mean(Psi_tilde)
  
  # influences
  Psi <- Psi_tilde - ate
  
  var <- mean(Psi^2)
  se <- sqrt(var / n)
  
  out <- c(table(D)[[2]], table(D)[[1]], ate, se)
  
  return(out)
}
