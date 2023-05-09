#' Auto DML based on rlasso
#'
#' Implements AutoDML based onlasso learners (\code{rlasso}) as developed
#' in Chernozhukov et al. (2022) and Klosin and Vilgalys (2022).
#' Estimates the average treatment effect (ATE)
#' of a binary treatment variable, \eqn{D} on an outcome variable \eqn{Y} in presence of confounding variables \eqn{X} or the average partial derivative (APD) of a continuous treatment variable, \eqn{D} on an outcome variable \eqn{Y}.
#' The algorithm is based on the double ML
#' cross-fitting algorithm. This method was first introduced in Chernozhukov et
#' al. (2018) which required manual calculation of the Riesz representer in the
#' corresponding (IRM) model. In AutoDML (Chernozhukov 2022), the Riesz
#' representer is estimated automatically and doesn't need to be explicitly
#' computed.
#'
#' 
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E.,
#' Hansen, C., Newey, W. and Robins, J. (2018a), Double/debiased machine learning
#' for treatment and structural parameters.
#' The Econometrics Journal, 21: C1-C68. \doi{10.1111/ectj.12097}.
#'
#' @references Chernozhukov, V., Newey, W.K., and Singh, R. (2022), Automatic
#' Debiased Machine Learning of Causal and Structural Effects.
#' Econometrica, 90.3: 967-1027
#' 
#' @references Klosin, S., Vilgalys, M. (2022), Estimating Continuous
#' Treatment Effects in Panel Data using Machine Learning with an Agricultural
#' Application. arXiv preprint arXiv:2207.08789 (2022).
#'
#' @param data Backend with preprocessed data for estimation of 
#' ATE (output of function DataATEAutoDML) or APD (output of DataAPDAutoDML function). 
#' @param debiased debiased (if TRUE) vs. biased (if FALSE) results.
#' Default is TRUE.
#' @param D_LB lower bound on D (default 0)
#' @param D_add value added to D (default 0.2)
#' @param L number of folds data is split into (default 5)
#' @param max_iter maximum iterations of Lasso for debiasing (default 10)
#' @param est_type specifying the type of estimator, must either be "APD" 
#' for average partial derivative or "ATE" for average treatment effect.
#' @param weights weight vector - note for panel must be length (N-1)*T - because weight vector
#' for first differenced data 
#' @return an object of class \code{rlassoAutoDML} with estimated effects,
#' standard errors and individual effects in the form of a \code{list}.
#' @seealso rlasso
#' @examples
#' library(hdm)
#' data <- DGP_contD(N, 20)
#' Xnames = colnames(data)[grep("X", colnames(data))]
#' backend_contD = DataAPDAutoDML(x = Xnames, d = "D", y = "y", data = data,
#'                               poly = 3, interactions = TRUE,
#'                               intercept = "intercept")
#' auto_apd = rlassoAutoDML(backend_contD, est_type = "APD", debiased = TRUE)
#' auto_apd$te
#' 
#' @export
#' @rdname rlassoAutoDML
rlassoAutoDML <- function(x, ...) {
  UseMethod("rlassoAutoDML")
} 

#' @export
#' @rdname rlassoAutoDML
rlassoAutoDML.default <- function(data, D_LB = 0, D_add = 0.2,
                                  debiased = TRUE, L = 5, max_iter = 10,
                                  gamma_learner = "rlasso",
                                  est_type = NULL, prelim_est = FALSE, 
                                  weights = NULL) {
  checkmate::check_subset(est_type, c("ATE", "APD"))
  if (est_type == "ATE") {
    checkmate::check_class(data, "DataATEAutoDML")
  }
  if (est_type == "APD") {
    checkmate::check_class(data, "DataAPDAutoDML")
  }
  checkmate::check_subset(gamma_learner, c("rlasso", "cv.glmnet"))
  checkmate::check_logical(prelim_est)

  if (prelim_est & est_type == "APD") {
    stop("Preliminary estimator only implemented for est_type = 'ATE'.")
  }
  
  p <- data$p
  n <- data$n
  Y <- data$Y
  M <- data$M
  X <- data$X
  
  if (prelim_est) {
    M_prelim <- data$M_prelim
    X_prelim <- data$X_prelim
  }
  prelim_quantities.nl <- NULL
  
  folds <- split(sample(n, n, replace = FALSE), as.factor(1:L))
  Psi_tilde <- rep(NA_real_, n)

  for (l in 1:L) {
    Y.l <- Y[folds[[l]]]
    Y.nl <- Y[-folds[[l]]]
    
    X.l <- X[folds[[l]], ]
    X.nl <- X[-folds[[l]], ]
    
    M.l <- M[folds[[l]], ]
    M.nl <- M[-folds[[l]], ]
    
    if (prelim_est) {
      M_prelim.nl <- M_prelim[-folds[[l]], ]
      X_prelim.nl <- X_prelim[-folds[[l]], ]
      prelim_quantities.nl <- list("M_prelim" = M_prelim.nl,
                                   "X_prelim" = X_prelim.nl,
                                   "p_prelim_out" = data$p_prelim_out)
    }
    
    # get stage 1 (on nl)
    # TODO: pass additional arguments via ...
    rho_hat <- RMD_stable(Y = Y.nl, X = X.nl, M = M.nl, p, D_LB, D_add,
                          max_iter, prelim_quantities = prelim_quantities.nl)
    
    if (gamma_learner == "rlasso") {
      gamma_coeff <- rlasso(X.nl, Y.nl, intercept = TRUE)$coefficients[-1]
        #rlasso(X.nl, Y.nl, intercept = FALSE)$coefficients
    } else {
      gamma_coeff <- coef(glmnet::cv.glmnet(x = X.nl, y = Y.nl, nfolds = 5,
                                            intercept = TRUE), lambda = "s.min")[-1]
    }
    
    print(paste0("fold: ", l))
    
    # get stage 2 (on l)
    # psi_star
    Psi_tilde.l <- psi_tilde(Y.l, X.l, M.l, rho_hat,
                             gamma_coeff, debiased)
    Psi_tilde[folds[[l]]] <- Psi_tilde.l
  }
  
  # point estimation
  if(is.null(weights)){
    ate <- mean(Psi_tilde)
  }else{
    ate <- weighted.mean(m_fit, weights)
    }
  
  
  # influences
  Psi <- Psi_tilde - ate
  
  var <- mean(Psi^2)
  se <- sqrt(var / n)
  
  res <- list(
    se = se,
    te = ate,
    individual = Psi_tilde,
    type = est_type,
    call = match.call(),
    samplesize = n)
  
  class(res) <- "rlassoAutoDML"
  
  return(res)
}


################# Methods for rlassoAutoDML

#' Methods for S3 object \code{rlassoAutoDML}
#'
#' Objects of class \code{rlassoAutoDML} are constructed by  \code{rlassoAutoDML}.
#' \code{print.rlassoAutoDML} prints and displays some information about fitted \code{rlassoAutoDML} objects.
#' \code{summary.rlassoAutoDML} summarizes information of a fitted \code{rlassoAutoDML} object.
#' \code{confint.rlassoAutoDML} extracts the confidence intervals.
#' @param object an object of class \code{rlassoAutoDML}
#' @param digits number of significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @rdname methods.rlassoAutoDML
#' @aliases methods.rlassoAutoDML print.rlassoAutoDML summary.rlassoAutoDML
#' @export

print.rlassoAutoDML <- function(object, digits = max(3L, getOption("digits") - 3L),
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

#' @rdname methods.rlassoAutoDML
#' @export
summary.rlassoAutoDML <- function(object, digits = max(3L, getOption("digits") -
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


#' @rdname methods.rlassoAutoDML
#' @export
confint.rlassoAutoDML <- function(object, parm, level = 0.95, ...) {
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
