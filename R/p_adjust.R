#'Multiple Testing Adjustment of p-values for S3 objects \code{rlassoEffects}
#'and \code{lm}
#'
#'Multiple hypotheses testing adjustment of p-values from a high-dimensional
#'linear model.
#'
#'Multiple testing adjustment is performed for S3 objects of class
#'\code{rlassoEffects} and \code{lm}. Implemented methods for multiple testing
#'adjustment are Romano-Wolf stepdown '\code{RW}' (default) and the adjustment
#'methods available in the \code{p.adjust} function of the \code{stats} package,
#'including the Bonferroni, Bonferroni-Holm, and Benjamini-Hochberg corrections,
#'see \code{\link{p.adjust.methods}}.
#'
#'Objects of class \code{rlassoEffects} are constructed by
#'\code{\link{rlassoEffects}}.
#'
#'@param x an object of S3 class \code{rlassoEffects} or \code{lm}.
#'@param method the method of p-value adjustment for multiple testing.
#'  Romano-Wolf stepdown ('\code{RW}') is chosen by default.
#'@param test.index vector of integers, logicals or variables names indicating
#'  the position of coefficients (integer case), logical vector of length of the
#'  coefficients (TRUE or FALSE) or the coefficient names of x which should be
#'  tested simultaneously (only for S3 class \code{lm}). If missing, all
#'  coefficients are considered.
#'@param B number of bootstrap repetitions (default 1000).
#'@param ... further arguments passed on to methods.
#'@rdname p_adjust
#'@aliases p_adjust.rlassoEffects p_adjust.lm
#'@return A matrix with the estimated coefficients and the p-values that are
#'  adjusted according to the specified method.
#'@references J.P. Romano, M. Wolf (2005). Exact and approximate stepdown
#'  methods for multiple hypothesis testing. Journal of the American Statistical
#'  Association, 100(469), 94-108.
#'@references J.P. Romano, M. Wolf (2016). Efficient computation of adjusted
#'  p-values for resampling-based stepdown multiple testing. Statistics and
#'  Probability Letters, (113), 38-40.
#'@references A. Belloni, V. Chernozhukov, K. Kato (2015). Uniform
#'  post-selection inference for least absolute deviation regression and other
#'  Z-estimation problems. Biometrika, 102(1), 77-94.
#'
#'
#' @examples
#' library(hdm);
#' set.seed(1)
#' n = 100 #sample size
#' p = 25 # number of variables
#' s = 3 # nubmer of non-zero variables
#' X = matrix(rnorm(n*p), ncol=p)
#' colnames(X) <- paste("X", 1:p, sep="")
#' beta = c(rep(3,s), rep(0,p-s))
#' y = 1 + X%*%beta + rnorm(n)
#' data = data.frame(cbind(y,X))
#' colnames(data)[1] <- "y"
#' lasso.effect = rlassoEffects(X, y, index=c(1:20))
#' pvals.lasso.effect = p_adjust(lasso.effect, method = "RW", B = 1000)
#' ols = lm(y ~ -1 + X, data)
#' pvals.ols = p_adjust(ols, method = "RW", B = 1000)
#' pvals.ols = p_adjust(ols, method = "RW", B = 1000, test.index = c(1,2,5))
#' pvals.ols = p_adjust(ols, method = "RW", B = 1000, test.index = c(rep(TRUE, 5), rep(FALSE, p-5)))
#' @export

p_adjust = function(x, ...){
  UseMethod("p_adjust")
}


#' @describeIn p_adjust \code{\link{rlassoEffects}}.
#' @export
#'
p_adjust.rlassoEffects <- function(x, method = "RW", B = 1000, ...) {
  checkmate::checkClass(x, "rlassoEffects")
  checkmate::checkChoice(method, c("RW", stats::p.adjust.methods))

  n <- x$samplesize
  k <- length(x$coefficients)
  cf <- coef(x)

  pinit <- corr.padj <- pval <- vector(mode = "numeric", length = k)

  if (is.element(method, stats::p.adjust.methods)) {
    pval <- stats::p.adjust(x$pval, method = method, n = k)
  }

  if (method == "RW") {
    e <- x$residuals$e
    v <- x$residuals$v
    ev <- e * v
    Ev2 <- colMeans(v^2)
    Omegahat <- matrix(NA, ncol = k, nrow = k)
    for (j in 1:k) {
      for (l in 1:k) {
        Omegahat[j, l] = Omegahat[l, j] = 1/(Ev2[j] * Ev2[l]) * mean(ev[, j] * ev[, l])
      }
    }
    se <- sqrt(diag(Omegahat))

    Beta_i <- matrix(NA, ncol = k, nrow = B)
    for (i in 1:B) {
      Beta_i[i, ] <- MASS::mvrnorm(mu = rep(0, k), Sigma = Omegahat/n)
    }

    tstats <- cf/se
    stepdown.index <- order(abs(tstats), decreasing = TRUE)
    ro <- order(stepdown.index)

    for (s in 1:k) {
      if (s == 1) {
        sim <- apply(Beta_i, 1, function(z) max(abs(z)/se))
        pinit[s] <- pmin(1, (sum(sim >= abs(tstats[stepdown.index][s])))/B)
      }
      if (s > 1) {
        sim <- apply(Beta_i[, -stepdown.index[1:(s - 1)], drop = F], 1, function(z) max(abs(z)/se[-stepdown.index[1:(s - 1)]]))
        pinit[s] <- pmin(1, (sum(sim >= abs(tstats[stepdown.index][s])))/B)
      }
    }

      for (j in 1:k) {
        if (j == 1) {
          corr.padj[j] <- pinit[j]
        }

        if (j > 1) {
          corr.padj[j] <- max(pinit[j], corr.padj[j - 1])
        }
      }
      pval <- corr.padj[ro]
  }

  res <- as.matrix(cbind(cf, pval))
  colnames(res) <- c("Estimate.", "pval")

  return(res)
}


#' @describeIn p_adjust \code{\link[stats]{lm}}.
#' @export
p_adjust.lm = function(x, method = "RW", B = 1000, test.index = NULL, ...) {
  checkmate::checkClass(x, "lm")
  checkmate::checkChoice(method, c("RW", stats::p.adjust.methods))
  checkmate::assert(checkmate::checkNull(test.index), checkmate::checkLogical(test.index), checkmate::checkNumeric(test.index), checkmate::checkCharacter(test.index))

  cf = coef(x)
  pnames = names(cf)

  if (is.null(test.index)) {
    k <- length(coef(x))
    index <- c(1:k)
  } else {
    if (is.logical(test.index)) {
      k <- sum(test.index)
      stopifnot(length(test.index) == length(coef(x)))
      index <- which(test.index == T)
    } else {

      if (is.numeric(test.index)) {
        index <- as.integer(test.index)
        stopifnot(all(test.index <= length(coef(x))) && length(test.index) <= length(coef(x)))
        k <- length(test.index)
      }

      if (is.character(test.index)) {
        stopifnot(all(is.element(test.index, pnames)))
        index <- which(is.element(pnames, test.index))
        k <- length(test.index)
      }
    }
  }

  if (is.element(method, stats::p.adjust.methods)) {
    pval <- stats::p.adjust(summary(x)$coefficients[index, "Pr(>|t|)"], method = method, n = k)
  }

  if (method == "RW") {
    tstats <- summary(x)$coefficients[index, c("t value")]
    pinit <- corr.padj <- pval <- vector(mode = "numeric", length = k)
    Omegahat <- vcov(x)[index, index]
    se <- sqrt(diag(Omegahat))
    Beta_i <- matrix(NA, ncol = k, nrow = B)

    for (i in 1:B) {
      Beta_i[i, ] <- MASS::mvrnorm(mu = rep(0, k), Sigma = Omegahat)
    }

    stepdown.index <- order(abs(tstats), decreasing = TRUE)
    ro <- order(stepdown.index)

    for (s in 1:k) {
      if (s == 1) {
        sim <- apply(Beta_i, 1, function(z) max(abs(z)/se))
        pinit[s] <- pmin(1,  (sum(sim >= abs(tstats[stepdown.index][s])))/B)
      }
      if (s > 1 ) {
        sim <- apply(Beta_i[, -stepdown.index[1:(s - 1)], drop = FALSE], 1, function(z) max(abs(z)/se[-stepdown.index[1:(s - 1)]]))
        pinit[s] <- pmin(1, (sum(sim >= abs(tstats[stepdown.index][s])))/B)
      }
    }

    for (j in 1:k) {
      if (j == 1) {
        corr.padj[j] <- pinit[j]
      }

      if (j > 1) {
        corr.padj[j] <- max(pinit[j], corr.padj[j - 1])
      }
    }

    pval <- corr.padj[ro]
  }

  res <- as.matrix(cbind(cf[index], pval))
  colnames(res) <- c("Estimate.", "pval")

  return(res)
}

