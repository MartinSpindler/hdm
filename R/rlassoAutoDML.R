#' Auto DML based on rlasso
#' 
#' TODO:
#' insert a meaningful description and include references from the literature.
#'
#' @param Y A vector of outputs
#' @param D A vector of treatment values
#' @param X A matrix of covariates
#' @param dict A dictionary
#' TODO: explain what the dictionary does
#' @param bias debiased vs. biased results
#' TODO: do users really need the biased results or are they just used for
#' illustration purposes
#' @param D_LB 
#' TODO: documentation for this parameter
#' @param D_add 
#' TODO: documentation for this parameter
#' @param L
#' TODO: documentation for this parameter
#' @param max_iter
#' TODO: documentation for this parameter
#' @return list with average treatment effect and standard error
#' @examples
#' rlassoDML(Y, T, X, dict)
#' @export
#' @rdname rlassoDML
rlassoAutoDML <- function(Y, D, X, dict = NULL, D_LB = 0, D_add = 0.2,
                          bias = FALSE, L = 5, max_iter = 10) {
  
  if (is.null(dict)) {
    # TODO: insert default behavior if no dictionary is provided by the user
    dict = function(x, y) {
      return(x)
    }
    
  }
  p <- length(dict(D[1], X[1, ]))

  # p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  # TODO: is this comment true or is p0 = dim(X0)/4 used?
  p0 <- ceiling(p / 4)
  if (p > 60) {
    p0 <- ceiling(p / 40)

  }
  n <- nrow(X)
  folds <- split(sample(n, n, replace = FALSE), as.factor(1:L))

  Psi_tilde <- numeric(0)

  for (l in 1:L) {

    # folds[[1]] = test_inds_global

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
    # print(n)
    # print(p)
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
