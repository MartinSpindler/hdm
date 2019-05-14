#' Instrumental Variable Estimation with Selection on the exogenous Variables by Lasso   
#'
#'
#' This function estimates the coefficient of an endogenous variable by employing Instrument Variables in a setting where the exogenous variables are high-dimensional and hence
#' selection on the exogenous variables is required.
#' The function returns an element of class \code{rlassoIVselectX}
#'
#' The implementation is a special case of of Chernozhukov et al. (2015).
#' The option \code{post=TRUE} conducts post-lasso estimation for the Lasso estimations, i.e. a refit of the
#' model with the selected variables. Exogenous variables 
#' \code{x} are automatically used as instruments and added to the
#' instrument set \code{z}.
#'
#' @param x exogenous variables in the structural equation (matrix)
#' @param d endogenous variables in the structural equation (vector or matrix)
#' @param y outcome or dependent variable in the structural equation (vector or matrix)
#' @param z set of potential instruments for the endogenous variables.
#' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
#' @param \dots arguments passed to the function \code{rlasso}
#' @return An object of class \code{rlassoIVselectX} containing at least the following
#' components: \item{coefficients}{estimated parameter vector}
#' \item{vcov}{variance-covariance matrix} \item{residuals}{
#' residuals} \item{samplesize}{sample size}
#' @references Chernozhukov, V., Hansen, C. and M. Spindler (2015). Post-Selection and Post-Regularization Inference in Linear
#' Models with Many Controls and Instruments
#' \emph{American Economic Review, Papers and Proceedings} 105(5), 486--490.
#' @export
#' @rdname rlassoIVselectX
#' @examples
#' library(hdm)
#' data(AJR); y = AJR$GDP; d = AJR$Exprop; z = AJR$logMort
#' x = model.matrix(~ -1 + (Latitude + Latitude2 + Africa + 
#'                            Asia + Namer + Samer)^2, data=AJR)
#' dim(x)
#'   #AJR.Xselect = rlassoIV(x=x, d=d, y=y, z=z, select.X=TRUE, select.Z=FALSE)
#'   AJR.Xselect = rlassoIV(GDP ~ Exprop +  (Latitude + Latitude2 + Africa + Asia + Namer + Samer)^2 |
#'              logMort +  (Latitude + Latitude2 + Africa + Asia + Namer + Samer)^2,
#'              data=AJR, select.X=TRUE, select.Z=FALSE)
#' summary(AJR.Xselect)
#' confint(AJR.Xselect)
rlassoIVselectX <- function(x, ...)
  UseMethod("rlassoIVselectX") # definition generic function

#' @export
#' @rdname rlassoIVselectX
rlassoIVselectX.default <- function(x,d,y,z, post=TRUE, ...) {
  d <- as.matrix(d)
  z <- as.matrix(z)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:ncol(d), sep="")
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:ncol(x), sep="")
  if (is.null(colnames(z)) & !is.null(z)) colnames(z) <- paste("z", 1:ncol(z), sep="")
  n <- length(y)
  numIV <- dim(z)[2]
  Z <- cbind(z,x)
  lasso.d.x <- rlasso(d ~ x, post=post, ...)
  Dr <- d - predict(lasso.d.x)
  lasso.y.x <- rlasso(y ~ x, post=post, ...)
  Yr <- y - predict(lasso.y.x)
  Zr <- matrix(NA, nrow=n, ncol=numIV)
  
  k <- 1 + dim(d)[2] + numIV 
  selection.matrix <- matrix(NA, ncol = k, nrow = dim(x)[2])
  colnames(selection.matrix) <- c("y", colnames(d),  colnames(z))
  rownames(selection.matrix) <- colnames(x)
  selection.matrix[,1] <- lasso.y.x$index
  selection.matrix[,2] <- lasso.d.x$index
  
  for (i in seq(length.out=numIV)) {
  lasso.z.x <- rlasso(z[,i] ~ x, post=post, ...)
  Zr[,i] <- z[,i] - predict(lasso.z.x)
  selection.matrix[, i + 2] <- lasso.z.x$index
  }
  result <- tsls(y = Yr,d = Dr,x=NULL, z = Zr, intercept=FALSE)
  coef <- as.vector(result$coefficient)
  se <- diag(sqrt(result$vcov))
  vcov <- result$vcov
  
  names(coef) <- names(se) <- colnames(d)
  res <- list(coefficients=coef, se=se, vcov=vcov, call=match.call(), samplesize=n, 
              selection.matrix = selection.matrix)
  class(res) <- "rlassoIVselectX"
  return(res)
}


#' @rdname rlassoIVselectX
#' @export
#' @param formula An object of class \code{Formula} of the form " y ~ x + d | x + z" with y the outcome variable,
#' d endogenous variable, z instrumental variables, and x exogenous variables.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), typically the environment from which \code{rlassoIVselectX} is called.
rlassoIVselectX.formula <- function(formula, data, post=TRUE, ...) {
  mat <- f.formula(formula, data, all.categories = FALSE)
  y <- mat$Y
  x <- mat$X
  d <- mat$D
  z <- mat$Z
  res <- rlassoIVselectX(x=x, d=d, y=y, z=z, post=post, ...)
  res$call <- match.call()
  return(res)
}


################# Methods for rlassoIVselectX

#' Methods for S3 object \code{rlassoIVselectX}
#'
#' Objects of class \code{rlassoIVselectX} are constructed by \code{rlassoIVselectX}. 
#' \code{print.rlassoIVselectX} prints and displays some information about fitted \code{rlassoIVselectX} objects.
#' \code{summary.rlassoIVselectX} summarizes information of a fitted \code{rlassoIVselectX} object.
#' \code{confint.rlassoIVselectX} extracts the confidence intervals.
#' @param object an object of class \code{rlassoIVselectX}
#' @param x an object of class \code{rlassoIVselectX}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level	the confidence level required.
#' @rdname methods.rlassoIVselectX
#' @aliases methods.rlassoIVselectX print.rlassoIVselectX summary.rlassoIVselectX
#' @export

print.rlassoIVselectX <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(coef(x))
}

#' @rdname methods.rlassoIVselectX
#' @export

summary.rlassoIVselectX <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficient)
    table <- matrix(NA,ncol=4,nrow=k)
    rownames(table) <- names(object$coefficient)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[,1] <- object$coefficient
    table[,2] <- sqrt(diag(as.matrix(object$vcov)))
    table[,3] <- table[,1]/table[,2]
    table[,4] <- 2*pnorm(-abs(table[,3]))
    print("Estimation and significance testing of the effect of target variables in the IV regression model")
    printCoefmat(table, digits=digits,  P.values=TRUE, has.Pvalue=TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoIVselectX
#' @export

confint.rlassoIVselectX <- function(object, parm, level=0.95, ...) {
  n <- object$samplesize
  k <- length(object$coefficients)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  #fac <- qt(a, n-k)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                             pct))
  #ses <- sqrt(diag(object$vcov))[parm]
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}

#' Coefficients from S3 objects \code{rlassoIVselectX}
#'
#' Method to extract coefficients and selection matrix from objects of class \code{rlassoIVselectX}.
#' 
#' Printing coefficients and selection matrix for S3 object \code{rlassoIVselectX}. The first column of the selection matrix reports the selection index for the lasso regression of \code{y} on \code{x} in the specified
#' \code{rlassoIVselectX} command. \code{"x"} indicates that a variable has been selected, i.e., the corresponding estimated coefficient is different from zero.
#' The second column contains the selection index for the lasso regression of \code{d} on \code{x} and the remaining columns
#' the index of selected variables \code{x} for the instruments \code{z}. The very last column collects all variables that have been selected in at least one of the lasso regressions. 
#' 
#' @param object an object of class \code{rlassoIVselectX}, usually a result of a call 
#' \code{rlassoIVselectX} or \code{rlassoIV} with options \code{select.X=TRUE} and
#' \code{select.Z=FALSE}.
#' @param selection.matrix if TRUE, a selection matrix is returned that indicates the selected variables from each regression.
#' Default is set to FALSE. See section on details for more information. 
#' @param complete general option of the function \code{coef}.
#' @param ... further arguments passed to functions coef. 
#' @export
#' @rdname coef.rlassoIVselectX
#' @examples
#' library(hdm)
#' data(AJR); y = AJR$GDP; d = AJR$Exprop; z = AJR$logMort
#' x = model.matrix(~ -1 + (Latitude + Latitude2 + Africa + 
#'                            Asia + Namer + Samer)^2, data=AJR)
#AJR.Xselect = rlassoIV(x=x, d=d, y=y, z=z, select.X=TRUE, select.Z=FALSE)
#' AJR.Xselect = rlassoIV(GDP ~ Exprop +  (Latitude + Latitude2 + Africa + Asia + Namer + Samer)^2 |
#'                          logMort +  (Latitude + Latitude2 + Africa + Asia + Namer + Samer)^2,
#'                        data=AJR, select.X=TRUE, select.Z=FALSE)
#' coef(AJR.Xselect) # Default behavior
#' coef(AJR.Xselect, selection.matrix = TRUE) # print selection matrix
coef.rlassoIVselectX <-  function(object, complete = TRUE, selection.matrix = FALSE, ...){
  
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
