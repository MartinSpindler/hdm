two.norm <- function(x) {
  return(sqrt(x %*% x))
}

psi_tilde <- function(y, x, m, rho, gamma, debiased) {
  # stage 2
  # fit alpha
  alpha_fit <- x %*% rho
  # fit/predict gamma
  gamma_fit <- x %*% gamma
  # fit derivative of gamma
  m_fit <- m %*% gamma
  if (debiased) {
    return(m_fit + alpha_fit * (y - gamma_fit))
  } else if (!debiased) {
    return(m_fit)
  }
}

default_dict_ATE <- function(d, z) {
  return(c(1, d, z))
}

get_D <- function(Y, X, M, rho_hat, p) {
  n <- nrow(X)
  df <- matrix(0, p, n)
  
  for (i in 1:n){
    df[,i] = as.numeric(X[i,]*as.vector(rho_hat %*% X[i,]) - M[i,])
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
#' TODO: Finalize documentation
#' @inheritParams rlassoAutoDML
#' @param c parameter to tune lambda (default 0.5)
#' @param gamma parameter to tune lambda (default 0.1)
#' @param tol minimum improvement to continue looping (default 1e-6)
#'
#' @export
RMD_stable <- function(Y, X, M, p, D_LB = 0, D_add = 0.2,
                       max_iter = 10, c = 0.5,
                       gamma = 0.1, tol = 1e-6, prelim_quantities = NULL) {

  rho_hat <- rep(0, p)
  n <- nrow(X)
  
  if (!is.null(prelim_quantities)) {
    M_prelim <- prelim_quantities$M_prelim
    X_prelim <- prelim_quantities$X_prelim
    M_hat_prelim <- colMeans(M_prelim)
    G_hat_prelim <- t(X_prelim)%*%X_prelim/nrow(X_prelim)
    
    # initial estimate
    rho_hat_prelim <- solve(G_hat_prelim, M_hat_prelim)
    rho_hat[1:prelim_quantities$p_prelim_out] <- rho_hat_prelim
  }
  
  # moments
  M_hat <- colMeans(M)
  G_hat <- t(X) %*% X/n
  
  # penalty
  lambda <- c * qnorm(1 - gamma / (2 * p)) / sqrt(n) # snippet

  ###########
  # alpha_hat
  ###########
  diff_rho <- 1
  k <- 1
  l <- 0.1

  while (diff_rho > tol & k <= max_iter) {
    # previous values
    rho_hat_old <- rho_hat

    # normalization
    D_hat_rho <- get_D(Y, X, M, rho_hat_old, p)
    D_hat_rho <- pmax(D_LB, D_hat_rho)
    D_hat_rho <- D_hat_rho + D_add

    L <- c(l, rep(1, p - 1)) # dictionary is ordered (constant,...)
    lambda_vec <- lambda * L * D_hat_rho # v3: insert D here
    rho_hat <- LassoShooting.fit(G_hat, M_hat, lambda_vec, XX = -G_hat/2,
                                 Xy = -M_hat/2, beta.start = rep(0, p))$coefficients
    # difference
    diff_rho <- two.norm(rho_hat - rho_hat_old)
    k <- k + 1
  }

  return(rho_hat)
}

RMD_opt <- function(Y, X, M, p, c = 0.5,
                        gamma = 0.1, prelim_quantities = NULL) {
  n <- nrow(X)
  # moments
  M_hat <- colMeans(M)
  G_hat <- t(X) %*% X/n
  
  # penalty
  lambda <- c * qnorm(1 - gamma / (2 * p)) / sqrt(n) # snippet
  # optimization setup 
  rho <- Variable(p)
  D_guess <- rep(1, p)
  np_L1 <- function(M_hat, G_hat, rho){
    -2*t(rho)%*%M_hat + CVXR::quad_form(rho, G_hat) #PSDWRAP(G_hat)
  }
  problem <- Problem(Minimize(np_L1(M_hat, G_hat, rho) + 
                                2*lambda*CVXR::norm1(CVXR::multiply(D_guess,rho))))
  
  result <- CVXR::solve(problem)
  rho_hat <- as.vector(result$getValue(rho))
  return(rho_hat)
}

