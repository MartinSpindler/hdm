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

psi_tilde <- function(y, d, z, m, rho, gamma, dict, debiased) {
  gamma_dz <- dict(d, z) %*% gamma

  if (debiased) {
    alpha_dz <- dict(d, z) %*% rho
    return(m(y, d, z, gamma, dict) + alpha_dz * (y - gamma_dz))
  } else if (!debiased) {
    return(m(y, d, z, gamma, dict))
  }
}

default_dict <- function(d, z) {
  return(c(1, d, z))
}

get_MNG <- function(Y, D, X, b) {
  p <- length(b(D[1], X[1, ]))
  n.nl <- length(D)

  B <- matrix(0, n.nl, p)
  M <- matrix(0, p, n.nl)
  N <- matrix(0, p, n.nl)

  for (i in 1:n.nl) {
    B[i, ] <- b(D[i], X[i, ])
    M[, i] <- b(1, X[i, ]) - b(0, X[i, ])
    N[, i] <- m2(Y[i], D[i], X[i, ], b) # this is a more general formulation for N
  }

  M_hat <- rowMeans(M)
  N_hat <- rowMeans(N)
  G_hat <- t(B) %*% B / n.nl

  return(list(M_hat, N_hat, G_hat, B))
}


get_D <- function(Y, D, X, rho_hat, b) {
  n <- nrow(X)
  p <- length(b(D[1], X[1, ]))

  df <- matrix(0, p, n)
  for (i in 1:n) {
    df[, i] <- b(D[i], X[i, ]) * as.vector(rho_hat %*% b(D[i], X[i, ])) - (b(1, X[i, ]) - b(0, X[i, ]))
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
#' @inheritParams rlassoAutoDML
#' @param c parameter to tune lambda (default 0.5)
#' @param gamma parameter to tune lambda (default 0.1)
#' @param tol minimum improvement to continue looping (default 1e-6)
#'
#' @export
RMD_stable <- function(Y, D, X, p0, D_LB = 0, D_add = 0.2, max_iter = 10, b = NULL, c = 0.5, gamma = 0.1, tol = 1e-6) {

  if (is.null(b)) {
    b <- default_dict
  }

  k <- 1
  l <- 0.1

  p <- length(b(D[1], X[1, ]))
  n <- length(D)

  # low-dimensional moments
  # TODO: Check if order of columns in X really should play a role, in this version they do!
  # Alternative I: Random choice
  # col_indx <- sample(p, size = p0, replace = FALSE)
  # X0 <- X[, col_indx, drop = FALSE]
  # Alternative II: Based on preliminary screening as in rlasso
  X0 <- X[, 1:p0, drop = FALSE]
  MNG0 <- get_MNG(Y, D, X0, b)
  M_hat0 <- MNG0[[1]]
  N_hat0 <- MNG0[[2]]
  G_hat0 <- MNG0[[3]]

  # initial estimate
  rho_hat0 <- solve(G_hat0, M_hat0)
  rho_hat <- c(rho_hat0, rep(0, p - ncol(G_hat0)))
  beta_hat0 <- solve(G_hat0, N_hat0)
  beta_hat <- c(beta_hat0, rep(0, p - ncol(G_hat0)))

  # moments
  MNG <- get_MNG(Y, D, X, b)
  M_hat <- MNG[[1]]
  N_hat <- MNG[[2]]
  G_hat <- MNG[[3]]

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
    D_hat_rho <- get_D(Y, D, X, rho_hat_old, b)
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
