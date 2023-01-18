#' Auto DML based on rlasso
#'
#' Implements AutoDML based on rigorous lasso learner (\code{rlasso}) as developed
#' in Chernozhukov et al. (2018b). Estimates the average partial derivative of a continuous treatment variable, \eqn{D} on an outcome variable \eqn{Y} in presence of confounding variables \eqn{X}.The algorithm is based on the double ML
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
#' @param Y outcome variable / dependent variable
#' @param delta delta for partial derivative (amount of shift).
#' @param X.up X with upward shift for D
#' @param X.down X with downward shift for D
#' @param debiased debiased (if TRUE) vs. biased (if FALSE) results.
#' Default is TRUE.
#' @param D_LB lower bound on D (default 0)
#' @param D_add value added to D (default 0.2)
#' @param L number of folds data is split into (default 5)
#' @param max_iter maximum iterations of Lasso (default 10)
#' @return an object of class \code{rlassoAPDAutoDML} with estimated effects,
#' standard errors and individual effects in the form of a \code{list}.
#' @seealso rlasso
#' @examples
#' library(hdm)
# TODO: Add example
#' 
#' @export
#' @rdname rlassoAPDAutoDML
rlassoAPDAutoDML <- function(x, ...) {
  UseMethod("rlassoAPDAutoDML")
} # definition of generic method

#' @export
#' @rdname rlassoAPDAutoDML
rlassoAPDAutoDML.default <- function(X, Y, delta = NULL, X.up = NULL, X.down = NULL, D_LB = 0, D_add = 0.2,
                                  debiased = TRUE, L = 5, max_iter = 10) {
  
  # TODO: Asserts for X, X.low, and X.up

  p <- length(X[1, ])
  
  n <- nrow(X)
  folds <- split(sample(n, n, replace = FALSE), as.factor(1:L))

  Psi_tilde <- numeric(0)
  
  for (l in 1:L) {
    Y.l <- Y[folds[[l]]]
    Y.nl <- Y[-folds[[l]]]

    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    X.up.l=X.up[folds[[l]],]
    X.up.nl=X.up[-folds[[l]],]
    
    X.down.l=X.down[folds[[l]],]
    X.down.nl=X.down[-folds[[l]],]
    
    n.l <- nrow(X.l)
    n.nl <- nrow(X.nl)
    
    # get stage 1 (on nl)
    # TODO: pass additional arguments via ...
    rho_hat <- RMD_stable(Y.nl, NULL, X.nl, X.up.nl, X.down.nl, p, delta, D_LB, D_add, max_iter)
    
    gamma_coeff <- rlasso(X.nl, Y.nl, intercept = F)$coefficients
    
    print(paste0("fold: ", l))
    
    # get stage 2 (on l)
    # psi_star
    Psi_tilde.l <- rep(0, n.l)
    for (i in 1:n.l) {
      Psi_tilde.l[i] <- psi_tilde(Y.l[i], NULL, X.l[i, ], NULL, rho_hat, gamma_coeff, NULL, debiased, X.up.l[i, ], X.down.l[i, ], delta)
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
    type = "PD",
    call = match.call(),
    samplesize = n)
  
  class(res) <- "rlassoAPDAutoDML"
  
  return(res)
}


################# Methods for rlassoAPDAutoDML

#' Methods for S3 object \code{rlassoAPDAutoDML}
#'
#' Objects of class \code{rlassoAPDAutoDML} are constructed by  \code{rlassoAPDAutoDML}.
#' \code{print.rlassoAPDAutoDML} prints and displays some information about fitted \code{rlassoAPDAutoDML} objects.
#' \code{summary.rlassoAPDAutoDML} summarizes information of a fitted \code{rlassoAPDAutoDML} object.
#' \code{confint.rlassoAPDAutoDML} extracts the confidence intervals.
#' @param object an object of class \code{rlassoAPDAutoDML}
#' @param digits number of significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @rdname methods.rlassoAPDAutoDML
#' @aliases methods.rlassoAPDAutoDML print.rlassoAPDAutoDML summary.rlassoAPDAutoDML
#' @export

print.rlassoAPDAutoDML <- function(object, digits = max(3L, getOption("digits") - 3L),
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

#' @rdname methods.rlassoAPDAutoDML
#' @export
summary.rlassoAPDAutoDML <- function(object, digits = max(3L, getOption("digits") -
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


#' @rdname methods.rlassoAPDAutoDML
#' @export
confint.rlassoAPDAutoDML <- function(object, parm, level = 0.95, ...) {
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
