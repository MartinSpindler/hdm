#' rigorous Lasso for Linear Models: Inference
#'
#' Estimation and inference of (low-dimensional) target coefficients in a high-dimensional linear model.
#'
#' The functions estimates (low-dimensional) target coefficients in a high-dimensional linear model.
#' An application is e.g. estimation of a treatment effect \eqn{\alpha_0} in a
#' setting of high-dimensional controls. The user can choose between the so-called post-double-selection method and partialling-out.
#' The idea of the double selection method is to select variables by Lasso regression of
#' the outcome variable on the control variables and the treatment variable on
#' the control variables. The final estimation is done by a regression of the
#' outcome on the treatment effect and the union of the selected variables in
#' the first two steps. In partialling-out first the effect of the regressors on the outcome and the treatment variable is taken out by Lasso and then a regression of the residuals is conducted. The resulting estimator for \eqn{\alpha_0} is normal
#' distributed which allows inference on the treatment effect. It presents a wrap function for \code{rlassoEffect} 
#' which does inference for a single variable.
#'
#' @param x matrix of regressor variables serving as controls and potential
#' treatments. For \code{rlassoEffect} it contains only controls, for \code{rlassoEffects} both controls and potential treatments. For \code{rlassoEffects} it must have at least two columns.
#' @param y outcome variable (vector or matrix)
#' @param index vector of integers, logicals or variables names indicating the position (column) of
#' variables (integer case), logical vector of length of the variables (TRUE or FALSE) or the variable names of \code{x} which should be used for inference / as treatment variables.
#' @param method method for inference, either 'partialling out' (default) or 'double selection'. 
#' @param I3 For the 'double selection'-method the logical vector \code{I3} has same length as the number of variables in \code{x};
#' indicates if variables (TRUE) should be included in any case to the model and they are exempt from selection. These variables should not be included in the \code{index}; hence the intersection with \code{index} must be the empty set.
#' In the case of partialling out it is ignored.
#' @param post logical, if post Lasso is conducted with default \code{TRUE}.
#' @param \dots parameters passed to the \code{rlasso} function.
#' @return The function returns an object of class \code{rlassoEffects} with the following entries: \item{coefficients}{vector with estimated
#' values of the coefficients for each selected variable} \item{se}{standard error (vector)}
#' \item{t}{t-statistic} \item{pval}{p-value} \item{samplesize}{sample size of the data set} \item{index}{index of the variables for which inference is performed}
#' @references A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on
#' treatment effects after selection among high-dimensional controls. The
#' Review of Economic Studies 81(2), 608-650.
#' @export
#' @rdname rlassoEffects
#' @examples
#' library(hdm); library(ggplot2)
#' set.seed(1)
#' n = 100 #sample size
#' p = 100 # number of variables
#' s = 3 # nubmer of non-zero variables
#' X = matrix(rnorm(n*p), ncol=p)
#' colnames(X) <- paste("X", 1:p, sep="")
#' beta = c(rep(3,s), rep(0,p-s))
#' y = 1 + X%*%beta + rnorm(n)
#' data = data.frame(cbind(y,X))
#' colnames(data)[1] <- "y"
#' fm = paste("y ~", paste(colnames(X), collapse="+"))
#' fm = as.formula(fm)                 
#' lasso.effect = rlassoEffects(X, y, index=c(1,2,3,50))
#' lasso.effect = rlassoEffects(fm, I = ~ X1 + X2 + X3 + X50, data=data)
#' print(lasso.effect)
#' summary(lasso.effect)
#' confint(lasso.effect)
#' plot(lasso.effect)
# library(hdm)
# ## DGP
# n <- 250
# p <- 100
# px <- 10
# X <- matrix(rnorm(n*p), ncol=p)
# beta <- c(rep(2,px), rep(0,p-px))
# intercept <- 1
# y <- intercept + X %*% beta + rnorm(n)
# ## fit rlassoEffects object with inference on three variables
# rlassoEffects.reg <- rlassoEffects(x=X, y=y, index=c(1,7,20))
# ## methods
# summary(rlassoEffects.reg)
# confint(rlassoEffects.reg, level=0.9)
rlassoEffects <- function(x, ...)
  UseMethod("rlassoEffects") # definition generic function 

#' @export
#' @rdname rlassoEffects
rlassoEffects.default <- function(x, y, index = c(1:ncol(x)), method = "partialling out", 
                                  I3 = NULL, post = TRUE, ...) {
  
  checkmate::checkChoice(method, c("partialling out", "double selection"))
  
  if (is.logical(index)) {
    k <- p1 <- sum(index)
  } else {
    k <- p1 <- length(index)
  }
  n <- dim(x)[1]
  x <- as.matrix(x)
  # preprocessing index numerischer Vektor
  if (is.numeric(index)) {
    index <- as.integer(index)
    stopifnot(all(index <= ncol(x)) && length(index) <= ncol(x))
  } else {
    # logical Vektor
    if (is.logical(index)) {
      stopifnot(length(index) == ncol(x))
      index <- which(index == T)
    } else {
      # character Vektor
      if (is.character(index)) {
        stopifnot(all(is.element(index, colnames(x))))
        index <- which(is.element(colnames(x), index))
      } else {
        stop("argument index has an invalid type")
      }
    }
  }
  if (method == "double selection") {
    # check validity of I3
    I3ind <- which(I3 == T)
    if (length(intersect(index, I3ind) != 0)) 
      stop("I3 and index must not overlap!")
  }
  
  if (is.null(colnames(x))) 
    colnames(x) <- paste("V", 1:dim(x)[2], sep = "")
  coefficients <- as.vector(rep(NA_real_, k))
  se <- rep(NA_real_, k)
  t <- rep(NA_real_, k)
  pval <- rep(NA_real_, k)
  lasso.regs <- vector("list", k)
  reside <- matrix(NA, nrow = n, ncol = p1)
  residv <- matrix(NA, nrow = n, ncol = p1)
  coef.mat <- list()
  selection.matrix <- matrix(NA, ncol = k, nrow = dim(x)[2])
  names(coefficients) <- names(se) <- names(t) <- names(pval) <- names(lasso.regs) <- colnames(reside) <- colnames(residv) <- colnames(selection.matrix) <- colnames(x)[index]
  rownames(selection.matrix) <- colnames(x)
  for (i in 1:k) {
    d <- x[, index[i], drop = FALSE]
    Xt <- x[, -index[i], drop = FALSE]
    I3m <- I3[-index[i]]
    lasso.regs[[i]] <- try(col <- rlassoEffect(Xt, y, d, method = method, 
                                               I3 = I3m, post = post, ...), silent = TRUE)
    if (class(lasso.regs[[i]]) == "try-error") {
      next
    } else {
      coefficients[i] <- col$alpha
      se[i] <- col$se
      t[i] <- col$t
      pval[i] <- col$pval
      reside[, i] <- col$residuals$epsilon
      residv[, i] <- col$residuals$v
      coef.mat[[i]] <- col$coefficients.reg
      selection.matrix[-index[i],i] <- col$selection.index
    }
  }
  #colnames(coef.mat) <- colnames(x)[index]
  residuals <- list(e = reside, v = residv)
  res <- list(coefficients = coefficients, se = se, t = t, pval = pval, 
              lasso.regs = lasso.regs, index = index, call = match.call(), samplesize = n, 
              residuals = residuals, coef.mat = coef.mat, selection.matrix = selection.matrix)
  class(res) <- "rlassoEffects"
  return(res)
}

#' @rdname rlassoEffects
#' @param formula An element of class \code{formula} specifying the linear model.
#' @param I An one-sided formula specifying the variables for which inference is conducted.
#' @param included One-sided formula of variables which should be included in any case (only for method="double selection").
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @export
rlassoEffects.formula <- function(formula, data, I, method = "partialling out", 
                                  included = NULL, post = TRUE, ...) {
  cl <- match.call()
  if (missing(data))  data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 1
  y <- model.response(mf, "numeric")
  n <- length(y)
  x <- model.matrix(mt, mf)[,-1, drop=FALSE]
  cn <- attr(mt, "term.labels")
  try(if (is.matrix(eval(parse(text=cn)))) cn <- colnames(eval(parse(text=cn))), silent=TRUE)
  I.c <- check_variables(I, cn)
  #I.c <- grep(cn[I.c],colnames(X))
  I.c <- which(colnames(x) %in% cn[I.c])
  I3 <- check_variables(included, cn)
  #I3 <- grep(cn[I.c],colnames(X))
  I3 <- which(colnames(x) %in% cn[I.c])
  #if (length(intersect(I.c, I3) != 0)) 
  #  stop("I and included should not contain the same variables!")
  
  est <- rlassoEffects(x, y, index = I.c, method = method, 
                       I3 = I3, post = post, ...)
  est$call <- cl
  return(est)
}

#' @rdname rlassoEffects
#' @param d variable for which inference is conducted (treatment variable)
#' @export
rlassoEffect <- function(x, y, d, method = "double selection", I3 = NULL, 
                         post = TRUE, ...) {
  d <- as.matrix(d, ncol = 1)
  y <- as.matrix(y, ncol = 1)
  kx <- dim(x)[2]
  n <- dim(x)[1]
  if (is.null(colnames(d))) 
    colnames(d) <- "d1"
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:kx, sep = "")
  if (method == "double selection") {
    I1 <- rlasso(d ~ x, post = post, ...)$index
    I2 <- rlasso(y ~ x, post = post, ...)$index
    
    
    if (is.logical(I3)) {
      I <- I1 + I2 + I3
      I <- as.logical(I)
    } else {
      I <- I1 + I2
      I <- as.logical(I)
      names(I) <- union(names(I1),names(I2))
    }
    if (sum(I) == 0) {
      I <- NULL
    }
    x <- cbind(d, x[, I, drop = FALSE])
    reg1 <- lm(y ~ x)
    alpha <- coef(reg1)[2]
    names(alpha) <- colnames(d)
    xi <- reg1$residuals * sqrt(n/(n - sum(I) - 1))
    if (is.null(I)) {
      reg2 <- lm(d ~ 1)
    }
    if (!is.null(I)) {
      reg2 <- lm(d ~ x[, -1, drop = FALSE])
    }
    v <- reg2$residuals
    var <- 1/n * 1/mean(v^2) * mean(v^2 * xi^2) * 1/mean(v^2)
    se <- sqrt(var)
    tval <- alpha/sqrt(var)
    pval <- 2 * pnorm(-abs(tval))
    if (is.null(I)) {
      no.selected <- 1
    } else {
      no.selected <- 0
    }
    res <- list(epsilon = xi, v = v)
    # results <- list(alpha=unname(alpha), se=drop(se), t=unname(tval),
    # pval=unname(pval), no.selected=no.selected,
    # coefficients=unname(alpha), coefficient=unname(alpha),
    # coefficients.reg=coef(reg1), residuals=res, call=match.call(),
    # samplesize=n)
    se <- drop(se)
    names(se) <- colnames(d)
    results <- list(alpha = alpha, se = se, t = tval, pval = pval, 
                    no.selected = no.selected, coefficients = alpha, coefficient = alpha, 
                    coefficients.reg = coef(reg1), selection.index = I, residuals = res, call = match.call(), 
                    samplesize = n)
  }
  
  if (method == "partialling out") {
    reg1 <- rlasso(y ~ x, post = post, ...)
    yr <- reg1$residuals
    reg2 <- rlasso(d ~ x, post = post, ...)
    dr <- reg2$residuals
    reg3 <- lm(yr ~ dr)
    alpha <- coef(reg3)[2]
    var <- vcov(reg3)[2, 2]
    se <- sqrt(var)
    tval <- alpha/sqrt(var)
    pval <- 2 * pnorm(-abs(tval))
    res <- list(epsilon = reg3$residuals, v = dr)
    I1 <- reg1$index
    I2 <- reg2$index
    I <- as.logical(I1 + I2)
    names(I) <- union(names(I1),names(I2))
    results <- list(alpha = unname(alpha), se = drop(se), t = unname(tval), 
                    pval = unname(pval), coefficients = unname(alpha), coefficient = unname(alpha), 
                    coefficients.reg = coef(reg1), selection.index = I, residuals = res, call = match.call(), 
                    samplesize = n)
  }
  class(results) <- "rlassoEffects"
  return(results)
}

################################################################################################################################### Methods for rlassoEffects

#' Methods for S3 object \code{rlassoEffects}
#'
#' Objects of class \code{rlassoEffects} are constructed by  \code{rlassoEffects}.
#' \code{print.rlassoEffects} prints and displays some information about fitted \code{rlassoEffect} objects.
#' summary.rlassoEffects summarizes information of a fitted \code{rlassoEffect} object and is described at \code{\link{summary.rlassoEffects}}.
#' \code{confint.rlassoEffects} extracts the confidence intervals.
#' \code{plot.rlassoEffects} plots the estimates with confidence intervals.
#'
#' @param x an object of class \code{rlassoEffects}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods.
#' @rdname methods.rlassoEffects
#' @aliases methods.rlassoEffects print.rlassoEffects confint.rlassoEffects plot.rlassoEffects
#' @export

print.rlassoEffects <- function(x, digits = max(3L, getOption("digits") - 
                                                  3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, 
                  quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(coef(x))
}

#' @rdname methods.rlassoEffects
#' @param object an object of class \code{rlassoEffects}
#' @param parm a specification of which parameters are to be given confidence intervals among the variables for which inference was done, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required
#' @param joint logical, if \code{TRUE} joint confidence intervals are calculated.
#' @export

# confint.rlassoEffects <- function(object, parm, level = 0.95, joint = FALSE, 
#                                   ...) {
#   B <- 500  # number of bootstrap repitions
#   n <- object$samplesize
#   k <- p1 <- length(object$coefficients)
#   cf <- coef(object)
#   pnames <- names(cf)
#   if (missing(parm)) 
#     parm <- pnames else if (is.numeric(parm)) 
#       parm <- pnames[parm]
#   if (!joint) {
#     a <- (1 - level)/2
#     a <- c(a, 1 - a)
#     # fac <- qt(a, n-k)
#     fac <- qnorm(a)
#     pct <- format.perc(a, 3)
#     ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
#                                                                pct))
#     ses <- object$se[parm]
#     ci[] <- cf[parm] + ses %o% fac
#   }
#   
#   if (joint) {
#     phi <- object$residuals$e * object$residuals$v
#     m <- 1/sqrt(colMeans(phi^2))
#     phi <- t(t(phi)/m)
#     sigma <- sqrt(colMeans(phi^2))
#     sim <- vector("numeric", length = B)
#     for (i in 1:B) {
#       xi <- rnorm(n)
#       phi_temp <- phi * xi
#       Nstar <- 1/sqrt(n) * colSums(phi_temp)
#       sim[i] <- max(abs(Nstar))
#     }
#     a <- (1 - level)/2
#     ab <- c(a, 1 - a)
#     pct <- format.perc(ab, 3)
#     ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
#                                                                pct))
#     hatc <- quantile(sim, probs = 1 - a)
#     ci[, 1] <- cf[parm] - hatc * 1/sqrt(n) * sigma
#     ci[, 2] <- cf[parm] + hatc * 1/sqrt(n) * sigma
#   }
#   return(ci)
# }


confint.rlassoEffects <- function(object, parm, level = 0.95, joint = FALSE, 
                                  ...) {
  B <- 500  # number of bootstrap repitions
  n <- object$samplesize
  k <- p1 <- length(object$coefficients)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames else if (is.numeric(parm)) 
      parm <- pnames[parm]
  if (!joint) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    # fac <- qt(a, n-k)
    fac <- qnorm(a)
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                               pct))
    ses <- object$se[parm]
    ci[] <- cf[parm] + ses %o% fac
  }
  
  if (joint) {
    e <- object$residuals$e
    v <- object$residuals$v
    ev <- e*v
    Ev2 <- colMeans(v^2)
    Omegahat <- matrix(NA, ncol=k, nrow=k)
    for (j in 1:k) {
      for (l in 1:k) {
        Omegahat[j,l] = Omegahat[l,j] = 1/(Ev2[j]*Ev2[l]) * mean(ev[,j]*ev[,l])
      }
    }
    var <- diag(Omegahat)
    names(var) <- names(cf)
    sim <- vector("numeric", length = B)
    for (i in 1:B) {
      beta_i <- MASS::mvrnorm(mu = rep(0,k), Sigma=Omegahat/n)
      sim[i] <- max(abs(beta_i/sqrt(var)))
    }
    a <- (1 - level) #not dividing by 2!
    ab <- c(a/2, 1 - a/2)
    pct <- format.perc(ab, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                               pct))
    hatc <- quantile(sim, probs = 1 - a)
    ci[, 1] <- cf[parm] - hatc * sqrt(var[parm])
    ci[, 2] <- cf[parm] + hatc * sqrt(var[parm])
  }
  return(ci)
}


#' @rdname methods.rlassoEffects
#' @export
#' @param main an overall title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param xlim vector of length two giving lower and upper bound of x axis
plot.rlassoEffects <- function(x, joint=FALSE, level= 0.95, main = "", xlab = "coef", ylab = "", 
                               xlim = NULL, ...) {
  
  # generate ordered KI-matrix
  coefmatrix <- cbind(summary(x)$coef, confint(x, joint = joint, level=level))[, c(1, 5, 6)]
  if (is.null(dim(coefmatrix))) {
    vec <- coefmatrix
    coefmatrix <- matrix(vec, ncol = 3)
    colnames(coefmatrix) <- names(vec)
  }
  
  rownames(coefmatrix) <- names(x$coefficients)
  coefmatrix <- as.data.frame(coefmatrix)
  coefmatrix <- cbind(rownames(coefmatrix), coefmatrix)
  colnames(coefmatrix) <- c("names", "coef", "lower", "upper")
  coefmatrix <- coefmatrix[order(abs(coefmatrix[, 2])), ]
  
  col <- "#000099"
  # scale
  if (missing(xlim)) {
    low <- min(coefmatrix[, -1])
    up <- max(coefmatrix[, -1])
  } else {
    low <- xlim[1]
    up <- xlim[2]
  }
  # generate points
  plotobject <- ggplot2::ggplot(coefmatrix, ggplot2::aes(y = coef, x = factor(names, 
                                                                              levels = names))) + ggplot2::geom_point(colour = col) + 
    ggplot2::geom_hline(colour = col, ggplot2::aes(width = 0.1, h = 0, yintercept=0))
  
  # generate errorbars (KIs)
  plotobject <- plotobject + ggplot2::geom_errorbar(ymin = coefmatrix$lower, 
                                                    ymax = coefmatrix$upper, colour = col)
  
  # further graphic parameter
  plotobject <- plotobject + ggplot2::ggtitle(main) + ggplot2::ylim(low, 
                                                                    up) + ggplot2::xlab(ylab) + ggplot2::ylab(xlab)
  
  
  ## invert x and y axis
  #plotobject <- plotobject + ggplot2::coord_flip()
  
  # layout
  plotobject <- plotobject + ggplot2::theme_bw() + ggplot2::geom_blank() + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())
  plotobject <- plotobject + ggplot2::scale_x_discrete(labels = abbreviate)
  # plot
  plotobject
}


################ Methods: summary

#' Summarizing rlassoEffects fits
#' 
#' Summary method for class \code{rlassoEffects}
#' 
#' Summary of objects of class \code{rlassoEffects}
#' 
#' @param object an object of class \code{rlassoEffects}, usually a result of a call to \code{rlassoEffects}
#' @param ... further arguments passed to or from other methods.
#' @rdname summary.rlassoEffects
#' @export 
summary.rlassoEffects <- function(object, ...) {
  ans <- NULL
  k <- length(object$coefficients)
  table <- matrix(NA, ncol = 4, nrow = k)
  rownames(table) <- names(object$coefficient)
  colnames(table) <- c("Estimate.", "Std. Error", "t value", "Pr(>|t|)")
  table[, 1] <- object$coefficients
  table[, 2] <- object$se
  table[, 3] <- object$t
  table[, 4] <- object$pval
  ans$coefficients <- table
  ans$object <- object
  class(ans) <- "summary.rlassoEffects"
  return(ans)
}


#' @param x an object of class \code{summary.rlassoEffects}, usually a result of a call or \code{summary.rlassoEffects}
#' @param digits the number of significant digits to use when printing.
#' @method print summary.rlassoEffects
#' @rdname summary.rlassoEffects
#' @export
print.summary.rlassoEffects <- function(x, digits = max(3L, getOption("digits") - 
                                                          3L), ...) {
  if (length(coef(x$object))) {
    k <- dim(x$coefficients)[1]
    table <- x$coefficients
    print("Estimates and significance testing of the effect of target variables")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}


#' Coefficients from S3 objects \code{rlassoEffects}
#'
#' Method to extract coefficients from objects of class \code{rlassoEffects}
#' 
#' Printing coefficients and selection matrix for S3 object \code{rlassoEffects}
#' 
#' @param object an object of class \code{rlassoEffects}, usually a result of a call \code{rlassoEffect} or \code{rlassoEffects}.
#' @param selection.matrix if TRUE, a selection matrix is returned that indicates the selected variables from each auxiliary regression. 
#' Default is set to FALSE. 
#' @param include.targets if FALSE (by default) only the selected control variables are listed in the \code{selection.matrix}. If set to TRUE, 
#' the selection matrix will also indicate the selection of the target coefficients that are specified in the  \code{rlassoEffects} call. 
#' @param complete general option of the function \code{coef}.
#' @param ... further arguments passed to functions coef or print. 
#' @rdname coef.rlassoEffects
#' @export
coef.rlassoEffects <- function(object, complete = TRUE, selection.matrix = FALSE, include.targets = FALSE, ...) {
  
  cf <- object$coefficients
  
  if (selection.matrix == TRUE) {
    
    mat <- object$selection.matrix
    
    if (is.null(mat)) {
      mat <- cbind(object$selection.index)
      dmat2 <- dim(mat)[2]
      rnames <- rownames(mat)
      targetindx <- stats::complete.cases(mat)
    }
    
    else {
      dmat2 <- dim(mat)[2]
      rnames <- rownames(mat)
      targetindx <- stats::complete.cases(mat)
      mat <- cbind(mat, as.logical(apply(mat, 1, sum)))
      colnames(mat)[dim(mat)[2]] <- "global"
    }
    
    if (include.targets == FALSE) {
      mat <- mat[targetindx, , drop = FALSE]
      rnames <- rownames(mat)
    }
    
    else{
      mat <- rbind(mat[targetindx == FALSE, , drop = FALSE], mat[targetindx, , drop = FALSE])
      rnames <- rownames(mat)
    }
    
    mat <- rbind(mat, apply(mat, 2, sum, na.rm = TRUE))
    mat <- apply(mat, 2, function(x) gsub(1, "x", x))
    mat <- apply(mat, 2, function(x) gsub(0, ".", x))
    mat[is.na(mat)] <- "-"
    rownames(mat) <- c(rnames, "sum")
    
    if (complete) {
      
      coef <- list(cf = cf, selection.matrix = mat)
      return(coef)
    }
    
    else {
      coef <- list(cf = cf[!is.na(cf)], selection.matrix = mat)
      return(coef)
    } 
  }
  
  else {
    
    if (complete) {
      return(cf)
    }
    
    else {
      return(cf[!is.na(cf)])
    } 
  }
  
}


#' Printing coefficients from S3 objects \code{rlassoEffects}
#'
#' Printing coefficients for class \code{rlassoEffects}
#' 
#' Printing coefficients and selection matrix for S3 object \code{rlassoEffects}
#' 
#' @param x an object of class \code{rlassoEffects}, usually a result of a call \code{rlassoEffect} or \code{rlassoEffects}.
#' @param selection.matrix if TRUE, a selection matrix is returned that indicates the selected variables from each auxiliary regression. 
#' Default is set to FALSE. 
#' @param include.targets if FALSE (by default) only the selected control variables are listed in the \code{selection.matrix}. If set to TRUE, 
#' the selection matrix will also indicate the selection of the target coefficients that are specified in the  \code{rlassoEffects} call. 
#' @param complete general option of the function \code{coef}.
#' @param ... further arguments passed to functions coef or print. 
#' @rdname print_coef
#' @aliases print_coef.rlassoEffects
#' @export
print_coef <-  function(x, ...){
  UseMethod("print_coef")
}


#' @rdname print_coef 
#' @export
print_coef.rlassoEffects <- function(x, complete = TRUE, selection.matrix = FALSE, include.targets = TRUE,  ...) {
  checkmate::check_class(x, "rlassoEffects")
  
  if (selection.matrix == FALSE) {
    cat("\n")
    print("Estimated target coefficients")
    print(coef(x), complete = complete, ...)
    cat("\n")
  }
  
  else {
    sel.mat <- coef(x, selection.matrix = selection.matrix, include.targets = include.targets, complete = complete, ...)
    cat("\n")
    print("Estimated target coefficients")
    print(sel.mat$cf)
    cat("\n")
    print("Matrix with selection index from auxiliary regressions")
    cat("\n")
    print(sel.mat$selection.matrix)
    cat("_ _ _ \n")
    print("'-' indicates a target variable; ")
    print("'x' indicates that a variable has been selected with rlassoEffects (coefficient is different from zero);") 
    print("'o' indicates that a variable has been de-selected with rlassoEffects (coefficient is zero).")
  }
}