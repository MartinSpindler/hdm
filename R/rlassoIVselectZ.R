#' Instrumental Variable Estimation with Lasso
#'
#' This function selects the instrumental variables in the first stage by
#' Lasso. First stage predictions are then used in the second stage as optimal
#' instruments to estimate the parameter vector. The function returns an element of class \code{rlassoIVselectZ}
#'
#' The implementation follows the procedure described in Belloni et al. (2012).
#' Option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the
#' model with the selected variables, to estimate the optimal instruments. The
#' parameter vector of the structural equation is then fitted by two-stage
#' least square (tsls) estimation.
#'
#' @param x exogenous variables in the structural equation (matrix)
#' @param d endogenous variables in the structural equation (vector or matrix)
#' @param y outcome or dependent variable in the structural equation (vector or matrix)
#' @param z set of potential instruments for the endogenous variables.
#' Exogenous variables serve as their own instruments.
#' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
#' @param intercept logical. If \code{TRUE}, intercept is included in the second stage equation.
#' @param \dots arguments passed to the function \code{rlasso}.
#' @return An object of class \code{rlassoIVselectZ} containing at least the following
#' components: \item{coefficients}{estimated parameter vector}
#' \item{vcov}{variance-covariance matrix} \item{residuals}{
#' residuals} \item{samplesize}{sample size} \item{selection.matrix}{matrix of selected variables in the first stage for each endogenous variable}
#' @references D. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
#' @export
rlassoIVselectZ <- function(x, ...)
  UseMethod("rlassoIVselectZ") # definition generic function

#' @export
#' @rdname rlassoIVselectZ
rlassoIVselectZ.default <- function(x, d, y, z, post = TRUE, intercept = TRUE, ...) {
  
  d <- as.matrix(d)
  if (is.vector(x)) x <- as.matrix(x)
  n <- length(y)
  kex <- dim(x)[2]
  ke <- dim(d)[2]
  
  if (is.null(colnames(d))) 
    colnames(d) <- paste("d", 1:ke, sep = "")
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:kex, sep = "")
  
  Z <- cbind(z,x) # including the x-variables as instruments
  kiv <- dim(Z)[2]
  select.mat <- NULL # matrix with the selected variables
  
  # first stage regression
  Dhat <- NULL
  flag.const <- 0
  for (i in 1:ke) {
    di <- d[, i]
    #lasso.fit <- rlasso(di ~ Z, post = post, intercept = intercept, ...)
    lasso.fit <- rlasso(y=di, x=Z, post = post, intercept = intercept, ...)
    if (sum(lasso.fit$ind) == 0) {
      dihat <- rep(mean(di), n)  #dihat <- mean(di)
      flag.const <- flag.const + 1
      if (flag.const >1) message("No variables selected for two or more instruments leading to multicollinearity problems.")
      #intercept <- FALSE # to avoid multicollineariry
      select.mat <- cbind(select.mat, FALSE)
    } else {
      # dihat <- z%*%lasso.fit$coefficients
      dihat <- predict(lasso.fit)
      select.mat <- cbind(select.mat, lasso.fit$index)
    }
    Dhat <- cbind(Dhat, dihat)
  }
  colnames(select.mat) <- colnames(d)
  #if (intercept) { #?
  #  Dhat <- cbind(Dhat, 1, x) 
  #  d <- cbind(d, 1, x)
  #} else {
  #  Dhat <- cbind(Dhat, x)
  #  d <- cbind(d, x)
  #}
  
  Dhat <- cbind(Dhat, x)
  d <- cbind(d, x)
  
  # calculation coefficients
  #alpha.hat <- solve(t(Dhat) %*% d) %*% (t(Dhat) %*% y)
  alpha.hat <- MASS::ginv(t(Dhat)%*%d)%*%(t(Dhat)%*%y)
  # calcualtion of the variance-covariance matrix
  residuals <- y - d %*% alpha.hat
  #Omega.hat <- t(Dhat) %*% diag(as.vector(residuals^2)) %*% Dhat  #  Dhat.e <- Dhat*as.vector(residuals);  Omega.hat <- t(Dhat.e)%*%Dhat.e
  Omega.hat <- t(Dhat) %*% (Dhat*as.vector(residuals^2))
  Q.hat.inv <- MASS::ginv(t(d) %*% Dhat)  #solve(t(d)%*%Dhat)
  vcov <- Q.hat.inv %*% Omega.hat %*% t(Q.hat.inv)
  rownames(alpha.hat) <- c(colnames(d))
  colnames(vcov) <- rownames(vcov) <- rownames(alpha.hat)
  
  if (is.null(x)){
  res <- list(coefficients = alpha.hat[1:ke, ], se = sqrt(diag(vcov))[1:ke], 
              vcov = vcov[1:ke, 1:ke, drop = FALSE], residuals = residuals, samplesize = n, selected = select.mat,
              selection.matrix =  select.mat,
              call = match.call()) 
  }
  else{
    if (is.null(x) == FALSE){
      res <- list(coefficients = alpha.hat[1:ke, ], se = sqrt(diag(vcov))[1:ke], 
                  vcov = vcov[1:ke, 1:ke, drop = FALSE], 
                  coefficients.controls = alpha.hat[(ke + 1):(ke + kex), ], se.controls = sqrt(diag(vcov))[(ke + 1):(ke + kex)], 
                  vcov.controls = vcov[(ke + 1):(ke + kex), (ke + 1):(ke + kex), drop = FALSE],
                  residuals = residuals, samplesize = n, selected = select.mat, selection.matrix =  select.mat,
                  call = match.call())
    }
  }
  class(res) <- "rlassoIVselectZ"
  return(res)
}


#' @rdname rlassoIVselectZ
#' @export
#' @param formula An object of class \code{Formula} of the form " y ~ x + d | x + z" with y the outcome variable,
#' d endogenous variable, z instrumental variables, and x exogenous variables.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), typically the environment from which \code{rlassoIVselectZ} is called.
rlassoIVselectZ.formula <- function(formula, data, post=TRUE, intercept = TRUE, ...) {
  mat <- f.formula(formula, data, all.categories = FALSE)
  y <- mat$Y
  x <- mat$X
  d <- mat$D
  z <- mat$Z
  
  res <- rlassoIVselectZ(x=x,d=d,y=y,z=z, post=post, intercept=intercept, ...)
  res$call <- match.call()
  return(res)
}
################# Methods for rlassoIVselectZ

#' Methods for S3 object \code{rlassoIVselectZ}
#'
#' Objects of class \code{rlassoIVselectZ} are constructed by \code{rlassoIVselectZ}. 
#' \code{print.rlassoIVselectZ} prints and displays some information about fitted \code{rlassoIVselectZ} objects.
#' \code{summary.rlassoIVselectZ} summarizes information of a fitted \code{rlassoIVselectZ} object.
#' \code{confint.rlassoIVselectZ} extracts the confidence intervals.
#' @param object an object of class \code{rlassoIVselectZ}
#' @param x an object of class \code{rlassoIVselectZ}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @rdname methods.rlassoIVselectZ
#' @aliases methods.rlassoIVselectZ print.rlassoIVselectZ summary.rlassoIVselectZ
#' @export

print.rlassoIVselectZ <- function(x, digits = max(3L, getOption("digits") - 
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

#' @rdname methods.rlassoIVselectZ
#' @export

summary.rlassoIVselectZ <- function(object, digits = max(3L, getOption("digits") - 
                                                           3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficients)
    table <- matrix(NA, ncol = 4, nrow = k)
    rownames(table) <- names(object$coefficients)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[, 1] <- object$coefficients
    table[, 2] <- sqrt(diag(as.matrix(object$vcov)))
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- 2 * pnorm(-abs(table[, 3]))
    print("Estimates and significance testing of the effect of target variables in the IV regression model")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoIVselectZ
#' @export

confint.rlassoIVselectZ <- function(object, parm, level = 0.95, ...) {
  n <- object$samplesize
  k <- length(object$coefficients)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames else if (is.numeric(parm)) 
      parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  # fac <- qt(a, n-k)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(as.matrix(object$vcov)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}


#' Coefficients from S3 objects \code{rlassoIVselectZ}
#'
#' Method to extract coefficients from objects of class \code{rlassoIVselectZ}.
#' 
#' Printing coefficients and selection matrix for S3 object \code{rlassoIVselectZ}. The columns of the selection matrix report the selection index for the first stage lasso regressions as specified
#' \code{rlassoIVselectZ} command, i.e., the selected variables for each of the endogenous variables. \code{"x"} indicates that a variable has been selected, i.e., the corresponding estimated coefficient is different from zero.
#' The very last column collects all variables that have been selected in at least one of the lasso regressions. 
#' 
#' @param object an object of class \code{rlassoIVselectZ}, usually a result of a call 
#' \code{rlassoIVselectZ} or \code{rlassoIV} with options \code{select.X=FALSE} and
#' \code{select.Z=TRUE}.
#' @param selection.matrix if TRUE, a selection matrix is returned that indicates the selected variables from each first stage regression.
#' Default is set to FALSE. See section on details for more information. 
#' @param complete general option of the function \code{coef}.
#' @param ... further arguments passed to functions coef. 
#' @export
#' @rdname coef.rlassoIVselectZ
#' @examples
#' \dontrun{
#' lasso.IV.Z = rlassoIVselectZ(x=x, d=d, y=y, z=z)
#' data(EminentDomain)
#' z <- EminentDomain$logGDP$z # instruments
#' x <- EminentDomain$logGDP$x # exogenous variables
#' y <- EminentDomain$logGDP$y # outcome varialbe
#' d <- EminentDomain$logGDP$d # treatment / endogenous variable
#' lasso.IV.Z = rlassoIVselectZ(x=x, d=d, y=y, z=z)
#' coef(lasso.IV.Z) # Default behavior
#' coef(lasso.IV.Z, selection.matrix = T)
#' }
coef.rlassoIVselectZ <-  function(object, complete = TRUE, selection.matrix = FALSE, ...){
  
  cf <- object$coefficients
  
  if (selection.matrix == TRUE) {
    
    mat <- object$selection.matrix
    
    dmat2 <- dim(mat)[2]
    rnames <- rownames(mat)
    mat <- cbind(mat, as.logical(apply(mat, 1, sum)))
    colnames(mat)[dim(mat)[2]] <- "global"
    mat <- rbind(mat, apply(mat, 2, sum, na.rm = TRUE))
    mat <- apply(mat, 2, function(x) gsub(1, "x", x))
    mat <- apply(mat, 2, function(x) gsub(0, ".", x))
    # mat[is.na(mat)] <- "-"
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


