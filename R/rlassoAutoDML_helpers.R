two.norm <- function(x) {
  return(sqrt(x %*% x))
}

m <- function(y, d, z, gamma, dict) { # all data arguments to make interchangeable with m2
  gamma_1z <- dict(1, z) %*% gamma
  gamma_0z <- dict(0, z) %*% gamma

  return(gamma_1z - gamma_0z)
}


m2 <- function(y, d, z, gamma) {
  return(y * gamma(d, z))
}

psi_tilde <- function(y, d, z, m, rho, gamma, dict, debiased, X.up = NULL, X.down = NULL, delta = NULL) {
  if (is.null(delta)) {
    gamma_dz <- dict(d, z) %*% gamma
    if (debiased) {
      alpha_dz <- dict(d, z) %*% rho
      return(m(y, d, z, gamma, dict) + alpha_dz * (y - gamma_dz))
    } else if (!debiased) {
      return(m(y, d, z, gamma, dict))
    }
  } else {
    gamma_dz <- dict(z, 0) %*% gamma
    gamma_dz_up <- dict(X.up, 0) %*% gamma
    gamma_dz_down <- dict(X.down, 0) %*% gamma
    if (debiased) {
      alpha_dz = dict(z, 0) %*% rho
      return((gamma_dz_up - gamma_dz_down)/delta + (alpha_dz - gamma_dz))
    } else if (!debiased) {
      return((gamma_dz_up - gamma_dz_down)/delta)
    }
  }
}


default_dict_ATE <- function(d, z) {
  return(c(1, d, z))
}

default_dict_PD <- function(d, z) {
  return(d)
}

get_MNG <- function(Y, D = NULL, X, b, p, X.up = NULL, X.down = NULL, delta = NULL) {
  
  n <- length(D)
  
  if (is.null(delta)) {
    if (is.null(p)) {
      p <- length(b(D[1], X[1, ]))
    }
    M <- matrix(0, p, n)
    N <- matrix(0, p, n)
    B <- matrix(0, n, p)
    
    for (i in 1:n) {
      B[i, ] <- b(D[i], X[i, ])
      M[, i] <- b(1, X[i, ]) - b(0, X[i, ])
      N[, i] <- m2(Y[i], D[i], X[i, ], b) # this is a more general formulation for N
    }
    M_hat <- rowMeans(M)
    N_hat <- rowMeans(N)
    G_hat <- t(B) %*% B / n
    
    return(list(M_hat = M_hat, N_hat = N_hat,
                G_hat = G_hat, B = B))
    
  } else {
    M <- matrix(0, p, n)
    N <- matrix(0, p, n)
    for (i in 1:n){
      M[,i]=(b(X.up[i,], X[i, ]) - b(X.down[i,], X[i,]))/delta #since m(w,b)=(x.up-x.down)/delta
      N[,i]=Y[i]*X[i,] #since m2(w,b)=y*x
    }
    M_hat <- rowMeans(M)
    N_hat <- rowMeans(N)
    G_hat = t(X) %*% X/n
    
    return(list(M_hat = M_hat, N_hat = N_hat, G_hat = G_hat))
  }
}


get_D <- function(Y, D = NULL, X, rho_hat, b, p, X.up = NULL, X.down = NULL, delta = NULL) {
  n <- nrow(X)
  df <- matrix(0, p, n)
  
  if (is.null(delta)) {
    
    for (i in 1:n) {
      df[, i] <- b(D[i], X[i, ]) * as.vector(rho_hat %*% b(D[i], X[i, ])) - (b(1, X[i, ]) - b(0, X[i, ]))
    }
  } else {
    
    for (i in 1:n){
      df[,i]=X[i,]*as.vector(rho_hat %*% X[i,])- (b(X.up[i,], X[i, ]) - b(X.down[i,], X[i, ]))/delta
    }
  }
  df <- df^2
  D2 <- rowMeans(df)
  D <- sqrt(D2)
  return(D) # pass around D as vector
}

#' Estimation of Riesz Representer
#'
#' Implements estimation of the Riesz representer \eqn{\alpha}.
#'
#' @inheritParams rlassoATEAutoDML
#' @param p0 initial dimension used in preliminary estimation step. By default, a rule of thumb is applied:
#' If \eqn{p \le 60}, p0 = p/4, else p0 = p/40.
#' @param c parameter to tune lambda (default 0.5)
#' @param gamma parameter to tune lambda (default 0.1)
#' @param tol minimum improvement to continue looping (default 1e-6)
#'
#' @export
RMD_stable <- function(Y, D = NULL, X = NULL, X.up = NULL, X.down = NULL, p, delta = NULL, D_LB = 0, D_add = 0.2, max_iter = 10, dict = NULL, p0 = NULL, c = 0.5, gamma = 0.1, tol = 1e-6) {

  if (is.null(p0)) {
    p0 <- ceiling(p / 4)
    if (p > 60) {
      p0 <- ceiling(p / 40)
    }
  }

  k <- 1
  l <- 0.1
  n <- length(D)

  # low-dimensional moments
  # TODO: Check if order of columns in X really should play a role, in this version they do!
  # Alternative I: Random choice
  # col_indx <- sample(p, size = p0, replace = FALSE)
  # X0 <- X[, col_indx, drop = FALSE]
  # Alternative II: Based on preliminary screening as in rlasso
  X0 <- X[, 1:p0, drop = FALSE]
  MNG0 <- get_MNG(Y, D, X0, dict, NULL, X.up, X.down, delta)
  M_hat0 <- MNG0$M_hat
  N_hat0 <- MNG0$N_hat
  G_hat0 <- MNG0$G_hat

  # initial estimate
  rho_hat0 <- solve(G_hat0, M_hat0)
  rho_hat <- c(rho_hat0, rep(0, p - ncol(G_hat0)))
  beta_hat0 <- solve(G_hat0, N_hat0)
  beta_hat <- c(beta_hat0, rep(0, p - ncol(G_hat0)))

  # moments
  MNG <- get_MNG(Y, D, X, dict, p, X.up, X.down, delta)
  M_hat <- MNG$M_hat
  N_hat <- MNG$N_hat
  G_hat <- MNG$G_hat

  # penalty
  lambda <- c * qnorm(1 - gamma / (2 * p)) / sqrt(n) # snippet

  ###########
  # alpha_hat
  ###########
  diff_rho <- 1

  while (diff_rho > tol & k <= max_iter) {
    # previous values
    rho_hat_old <- rho_hat + 0

    # normalization
    D_hat_rho <- get_D(Y, D, X, rho_hat_old, dict, p, X.up, X.down, delta)
    D_hat_rho <- pmax(D_LB, D_hat_rho)
    D_hat_rho <- D_hat_rho + D_add

    L <- c(l, rep(1, p - 1)) # dictionary is ordered (constant,...)
    lambda_vec <- lambda * L * D_hat_rho # v3: insert D here
    rho_hat <- LassoShooting.fit(G_hat, M_hat, lambda_vec, XX = -G_hat / 2, Xy = -M_hat / 2, beta.start = rep(0, p))$coefficients
    # difference
    diff_rho <- two.norm(rho_hat - rho_hat_old)
    k <- k + 1

  }

  return(rho_hat)
}
