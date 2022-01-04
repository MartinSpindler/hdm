#' Auto DML based on rlasso
#'
#' Implements AutoDML based on the Double ML algorithm for estimating causal effects.
#' This method was first introduced in Chernozhukov et al. (2018a) which required
#' manual calculation of the Riesz representer in the corresponding (IRM) model.
#' In AutoDML (Chernozhukov et al., 2018b), the Riesz representer is estimated
#' automatically and doesn't need to be explicitly computed. This implementation
#' is based on the \code{rlasso} learner.#'
#'
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E.,
#' Hansen, C., Newey, W. and Robins, J. (2018a), Double/debiased machine learning
#' for treatment and structural parameters.
#' The Econometrics Journal, 21: C1-C68. \doi{10.1111/ectj.12097}.
#'
#' @references Chernozhukov, V., Newey, W.K., and Singh, R. (2018b), Automatic
#' Debiased Machine Learning of Causal and Structural Effects. arXiv preprint
#' arXiv:1809.05224 (2018).
#'
#' @param Y outcome variable / dependent variable
#' @param D treatment variable (binary)
#' @param X exogenous variables
#' @param dict a dictionary
#' function of (d,z) that maps to a vector. Default is (1,d,z)
#' @param debiased debiased (if TRUE) vs. biased (if FALSE) results.
#' Default is TRUE.
#' @param D_LB lower bound on D (default 0)
#' @param D_add value added to D (default 0.2)
#' @param L number of folds data is split into (default 5)
#' @param max_iter maximum iterations of Lasso (default 10)
#' @param p0 initial value of p. By default, a rule of thumb is applied.
# TODO: Explicitly describe rule of thumb.
#' @return named list with average treatment effect and standard error
#' @examples
#' # TODO: Include replicable example
#' # data = simulate_data(500)
#' #
#' # Y = data[[1]]
#' # D = data[[2]]
#' # X = data[[3]]
#'
#' # rlassoAutoDML(Y,D,X,dict = b2)
#' # rlassoAutoDML(Y, T, X)
#' @export
#' @rdname rlassoAutoDML
rlassoAutoDML <- function(x, ...) {
  UseMethod("rlassoAutoDML")
} # definition of generic method

#' @export
#' @rdname rlassoAutoDML
rlassoAutoDML.default <- function(Y, D, X, dict = NULL, D_LB = 0, D_add = 0.2,
                                  debiased = TRUE, L = 5, max_iter = 10, p0 = NULL) {

  if (is.null(dict)) {
    dict <- default_dict
  }

  p <- length(dict(D[1], X[1, ]))

  if (is.null(p0)) {
    p0 <- ceiling(p / 4)
    if (p > 60) {
      p0 <- ceiling(p / 40)
    }
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

    alpha_hat <- dict_mult_coef(d, z, rho_hat)

    n <- nrow(X.nl)
    p <- length(dict(T.nl[1], X.nl[1, ]))
    # Apply the dictionary b to W
    B <- matrix(0, n, p)
    for (i in 1:n) {
      B[i, ] <- dict(T.nl[i], X.nl[i, ])
    }

    gamma_coeff <- rlasso(B, Y.nl, intercept = F)$coefficients
    gamma_hat <- dict_mult_coef(d, z, gamma_coef)

    print(paste0("fold: ", l))

    # get stage 2 (on l)
    # psi_star
    Psi_tilde.l <- rep(0, n.l)
    for (i in 1:n.l) {
      Psi_tilde.l[i] <- psi_tilde(Y.l[i], T.l[i], X.l[i, ], m, alpha_hat, gamma_hat, debiased)
    }
    Psi_tilde <- c(Psi_tilde, Psi_tilde.l)
  }

  # point estimation
  ate <- mean(Psi_tilde)

  # influences
  Psi <- Psi_tilde - ate

  var <- mean(Psi^2)
  se <- sqrt(var / n)

  res <- list(
    "treated:" = table(D)[[2]],
    "untreated" = table(D)[[1]],
    "ATE:" = ate,
    "SE:" = se
  )

  class(res) <- "rlassoAutoDML"

  return(res)
}


#' prints output of rlassoAutoDML in an easy to read format
#' @param obj Object of class rlassoAutoDML
#'
#' @export
print.rlassoAutoDML <- function(obj) {
  print(paste(" treated: ", obj$treated, " untreated: ", obj$untreated, "   ATE:    ", round(obj$ATE, 2), "   SE:   ", round(obj$SE, 2), sep = ""))
}
