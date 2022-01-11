#' Auto DML based on rlasso
#'
#' Implements AutoDML based on rigorous lasso learner (\code{rlasso}) as developed
#' in Chernozhukov et al. (2018b). Estimates the average treatment effect (ATE) of a binary treatment variable, \eqn{D} on an outcome variable \eqn{Y} in presence of confounding variables \eqn{X}.
#' The algorithm is based on the double ML
#' cross-fitting algorithm. This method was first introduced in Chernozhukov et
#' al. (2018a) which required manual calculation of the Riesz representer in the
#' corresponding (IRM) model. In AutoDML (Chernozhukov et al., 2018b), the Riesz
#' representer is estimated automatically and doesn't need to be explicitly
#' computed.
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
#' @param X exogenous variables
#' @param D treatment variable (binary)
#' @param Y outcome variable / dependent variable
#' @param dict a dictionary
#' function of (d,z) that maps to a vector. Default is (1,d,z)
#' @param debiased debiased (if TRUE) vs. biased (if FALSE) results.
#' Default is TRUE.
#' @param D_LB lower bound on D (default 0)
#' @param D_add value added to D (default 0.2)
#' @param L number of folds data is split into (default 5)
#' @param max_iter maximum iterations of Lasso (default 10)
#' @return an object of class \code{rlassoATEAutoDML} with estimated effects,
#' standard errors and individual effects in the form of a \code{list}.
#' @seealso rlasso
#' @examples
#' library(hdm)
#' library(DoubleML)
#' data <- make_irm_data(theta = 2, return_type = "matrix")
#' Y <- data$y
#' D <- data$d
#' X <- data$X
#' 
#' # Define a dictionary
#' b2 <- function(d,z){
#'   return(c(1,d,z,d*z))
#' }
#' lasso.auto <- rlassoATEAutoDML(X, D, Y, dict = b2)
#' summary(lasso.auto)
#' print(lasso.auto)
#' confint(lasso.auto)
#' 
#' @export
#' @rdname rlassoATEAutoDML
rlassoATEAutoDML <- function(x, ...) {
  UseMethod("rlassoATEAutoDML")
} # definition of generic method

#' @export
#' @rdname rlassoATEAutoDML
rlassoATEAutoDML.default <- function(X, D, Y, dict = NULL, D_LB = 0, D_add = 0.2,
                                  debiased = TRUE, L = 5, max_iter = 10) {
  
  # TODO: Asserts for D
  
  if (is.null(dict)) {
    dict <- default_dict_ATE
  }
  
  p <- length(dict(D[1], X[1, ]))
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
    rho_hat <- RMD_stable(Y.nl, T.nl, X.nl, NULL, NULL, p, NULL, D_LB, D_add, max_iter, dict)
    
    n <- nrow(X.nl)
    p <- length(dict(T.nl[1], X.nl[1, ]))
    # Apply the dictionary b to W
    B <- matrix(0, n, p)
    for (i in 1:n) {
      B[i, ] <- dict(T.nl[i], X.nl[i, ])
    }
    
    gamma_coeff <- rlasso(B, Y.nl, intercept = F)$coefficients
    
    print(paste0("fold: ", l))
    
    # get stage 2 (on l)
    # psi_star
    Psi_tilde.l <- rep(0, n.l)
    for (i in 1:n.l) {
      Psi_tilde.l[i] <- psi_tilde(Y.l[i], T.l[i], X.l[i, ], m, rho_hat, gamma_coeff, dict, debiased)
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
    se = se,
    te = ate,
    individual = Psi_tilde,
    type = "ATE",
    call = match.call(),
    samplesize = n,
    treated = table(D)[[2]],
    untreated = table(D)[[1]])
  
  class(res) <- "rlassoATEAutoDML"
  
  return(res)
}


################# Methods for rlassoATEAutoDML

#' Methods for S3 object \code{rlassoATEAutoDML}
#'
#' Objects of class \code{rlassoATEAutoDML} are constructed by  \code{rlassoATEAutoDML}.
#' \code{print.rlassoATEAutoDML} prints and displays some information about fitted \code{rlassoATEAutoDML} objects.
#' \code{summary.rlassoATEAutoDML} summarizes information of a fitted \code{rlassoATEAutoDML} object.
#' \code{confint.rlassoATEAutoDML} extracts the confidence intervals.
#' @param object an object of class \code{rlassoATEAutoDML}
#' @param digits number of significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @rdname methods.rlassoATEAutoDML
#' @aliases methods.rlassoATEAutoDML print.rlassoATEAutoDML summary.rlassoATEAutoDML
#' @export

print.rlassoATEAutoDML <- function(object, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
  )
  if (length(object$te)) {
    cat("Treatment Effect\n")
    cat(paste("Type:", object$type), "\n")
    cat("Value:\n")
    print.default(format(object$te, digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No treatment effect\n")
  }
  cat("\n")
  invisible(object$te)
}

#' @rdname methods.rlassoATEAutoDML
#' @export
summary.rlassoATEAutoDML <- function(object, digits = max(3L, getOption("digits") -
                                                         3L), ...) {
  if (length(object$te)) {
    table <- matrix(NA, ncol = 4, nrow = 1)
    rownames(table) <- "TE"
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[, 1] <- object$te
    table[, 2] <- object$se
    table[, 3] <- table[, 1] / table[, 2]
    table[, 4] <- 2 * pnorm(-abs(table[, 3]))
    cat(
      "Estimation and significance testing of the treatment effect",
      "\n"
    )
    cat(paste("Type:", object$type), "\n")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}


#' @rdname methods.rlassoATEAutoDML
#' @export
confint.rlassoATEAutoDML <- function(object, parm, level = 0.95, ...) {
  n <- object$samplesize
  k <- 1
  cf <- object$te
  pnames <- "ATE"
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  fac <- qt(a, n - k)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(pnames), 2L), dimnames = list(
    pnames,
    pct
  ))
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}
