#' Post-Selection and Post-Regularization Inference in Linear Models with Many
#' Controls and Instruments
#'
#' The function estimates a treatment effect in a setting with very many
#' controls and very many instruments (even larger than the sample size).
#'
#' The implementation for selection on x and z follows the procedure described in Chernozhukov et al.
#' (2015) and is built on 'triple selection' to achieve an orthogonal moment
#' function. The function returns an object of S3 class \code{rlassoIV}.
#' Moreover, it is wrap function for the case that selection should be done only with the instruments Z (\code{rlassoIVselectZ}) or with 
#' the control variables X (\code{rlassoIVselectX}) or without selection (\code{tsls}). Exogenous variables 
#' \code{x} are automatically used as instruments and added to the
#' instrument set \code{z}.
#'
#' @aliases rlassoIV rlassoIVmult
#' @param x matrix of exogenous variables
#' @param d endogenous variable
#' @param y outcome / dependent variable (vector or matrix)
#' @param z matrix of instrumental variables
#' @param post logical, wheter post-Lasso should be conducted (default=\code{TRUE})
#' @param select.Z logical, indicating selection on the instruments.
#' @param select.X logical, indicating selection on the exogenous variables.
#' @param \dots arguments passed to the function \code{rlasso}
#' @return an object of class \code{rlassoIV} containing at least the following
#' components: \item{coefficients}{estimated parameter value}
#' \item{se}{variance-covariance matrix}
#' @references V. Chernozhukov, C. Hansen, M. Spindler (2015). Post-selection
#' and post-regularization inference in linear models with many controls and
#' instruments. American Economic Review: Paper & Proceedings 105(5), 486--490.
#' @rdname rlassoIV
#' @export
#' @examples
#'\dontrun{
#' data(EminentDomain)
#' z <- EminentDomain$logGDP$z # instruments
#' x <- EminentDomain$logGDP$x # exogenous variables
#' y <- EminentDomain$logGDP$y # outcome varialbe
#' d <- EminentDomain$logGDP$d # treatment / endogenous variable
#' lasso.IV.Z = rlassoIV(x=x, d=d, y=y, z=z, select.X=FALSE, select.Z=TRUE) 
#' summary(lasso.IV.Z)
#' confint(lasso.IV.Z)
#' }
rlassoIV <- function(x, ...)
  UseMethod("rlassoIV") # definition generic function

#' @rdname rlassoIV
#' @export
rlassoIV.default <- function(x, d, y, z, select.Z = TRUE, select.X = TRUE, post = TRUE, 
                     ...) {
  d <- as.matrix(d)
  z <- as.matrix(z)
  if (is.null(colnames(d))) 
    colnames(d) <- paste("d", 1:ncol(d), sep = "")
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:ncol(x), sep = "")
  if (is.null(colnames(z)) & !is.null(z)) 
    colnames(z) <- paste("z", 1:ncol(z), sep = "")
  n <- length(y)
  
  if (select.Z == FALSE && select.X == FALSE) {
    res <- tsls(x, d, y, z, homoscedastic = FALSE, ...)
    return(res)
  }
  
  if (select.Z == TRUE && select.X == FALSE) {
    res <- rlassoIVselectZ(x, d, y, z, post = post, ...)
    return(res)
  }
  
  if (select.Z == FALSE && select.X == TRUE) {
    res <- rlassoIVselectX(x, d, y, z, post = post, ...)
    return(res)
  }
  
  if (select.Z == TRUE && select.X == TRUE) {
    
    Z <- cbind(z, x)
    lasso.d.zx <- rlasso(Z, d, post = post, ...)
    lasso.y.x <- rlasso(x, y, post = post, ...)
    lasso.d.x <- rlasso(x, d, post = post, ...)
    if (sum(lasso.d.zx$index) == 0) {
      message("No variables in the Lasso regression of d on z and x selected")
      return(list(alpha = NA, se = NA))
    }
    selection.matrixZ <- matrix(NA, ncol = dim(d)[2], nrow = dim(Z)[2])
    rownames(selection.matrixZ) <- colnames(Z)
    colnames(selection.matrixZ) <- colnames(d)
    selection.matrixZ[,1] <- ind.dzx <- lasso.d.zx$index
    
    #PZ <- Z[, ind.dzx] %*% MASS::ginv(t(Z[, ind.dzx]) %*% Z[, ind.dzx]) %*% 
    #  t(Z[, ind.dzx]) %*% d
    PZ <- as.matrix(predict(lasso.d.zx))
    lasso.PZ.x <- rlasso(x, PZ, post = post, ...)
    
    selection.matrix <- matrix(NA, ncol = (1 + dim(d)[2]), nrow = dim(x)[2])
    rownames(selection.matrix) <- colnames(x)
    colnames(selection.matrix) <- c("y", colnames(d))
    selection.matrix[ , 1] <- lasso.y.x$index
    selection.matrix[ , 2] <- ind.PZx <- lasso.PZ.x$index
    
    if (sum(ind.PZx) == 0) {
      Dr <- d - mean(d)
    } else {
      Dr <- d - predict(lasso.PZ.x)  #x[,ind.PZx]%*%MASS::ginv(t(x[,ind.PZx])%*%x[,ind.PZx])%*%t(x[,ind.PZx])%*%PZ
    }
    
    if (sum(lasso.y.x$index) == 0) {
      Yr <- y - mean(y)
    } else {
      Yr <- lasso.y.x$residuals
    }
    
    if (sum(lasso.PZ.x$index) == 0) {
      Zr <- PZ - mean(x)
    } else {
      Zr <- lasso.PZ.x$residuals
    }
    
    result <- tsls(y = Yr, d = Dr, x = NULL, z = Zr, intercept = FALSE, homoscedastic = FALSE)
    coef <- as.vector(result$coefficient)
    se <- diag(sqrt(result$vcov))
    names(coef) <- names(se) <- colnames(d)
    res <- list(coefficients = coef, se = se, vcov = vcov, call = match.call(), 
                samplesize = n, selection.matrixZ = selection.matrixZ, selection.matrix = selection.matrix)
    class(res) <- "rlassoIV"
    return(res)
  }
}


#' @param formula An object of class \code{Formula} of the form " y ~ x + d | x + z" with y the outcome variable,
#' d endogenous variable, z instrumental variables, and x exogenous variables.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), typically the environment from which \code{rlassoIV} is called.
#' @rdname rlassoIV
#' @export
rlassoIV.formula <- function(formula, data, select.Z = TRUE, select.X = TRUE, post = TRUE, 
                     ...) {
  
  mat <- f.formula(formula, data, all.categories = FALSE)
 
  y <- mat$Y
  x <- mat$X
  d <- mat$D
  z <- mat$Z
  
  res <- rlassoIV(x=x, d=d, y=y, z=z, select.Z = select.Z, select.X = select.X, post = post, 
                  ...)
  res$call <- match.call()
  return(res)
  
}


################# Methods for rlassoIV

#' Methods for S3 object \code{rlassoIV}
#'
#' Objects of class \code{rlassoIV} are constructed by \code{rlassoIV}. 
#' \code{print.rlassoIV} prints and displays some information about fitted \code{rlassoIV} objects.
#' \code{summary.rlassoIV} summarizes information of a fitted \code{rlassoIV} object.
#' \code{confint.rlassoIV} extracts the confidence intervals.
#' @param object An object of class \code{rlassoIV}
#' @param x an object of class \code{rlassoIV}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @rdname methods.rlassoIV
#' @aliases methods.rlassoIV print.rlassoIV summary.rlassoIV
#' @export

print.rlassoIV <- function(x, digits = max(3L, getOption("digits") - 3L), 
                           ...) {
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

#' @rdname methods.rlassoIV
#' @export

summary.rlassoIV <- function(object, digits = max(3L, getOption("digits") - 
                                                    3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficient)
    table <- matrix(NA, ncol = 4, nrow = k)
    rownames(table) <- names(object$coefficients)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[, 1] <- object$coefficients
    table[, 2] <- object$se
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- 2 * pnorm(-abs(table[, 3]))
    cat("Estimates and Significance Testing of the effect of target variables in the IV regression model", 
        "\n")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoIV
#' @export

confint.rlassoIV <- function(object, parm, level = 0.95, ...) {
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
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}

#' Coefficients from S3 objects \code{rlassoIV}
#'
#' Method to extract coefficients from objects of class \code{rlassoIV}.
#' 
#' Printing coefficients and selection matrix for S3 object \code{rlassoIV}. \code{"x"} indicates that a variable has been selected, i.e., the corresponding estimated coefficient is different from zero.
#' The very last column collects all variables that have been selected in at least one of the lasso regressions represented in the \code{selection.matrix}. 
#' \code{rlassoIV} performs three lasso regression steps. A first stage lasso regression of the endogenous treatment variable \code{d} on the instruments \code{z} and exogenous covariates \code{x},
#' a lasso regression of \code{y} on the exogenous variables \code{x}, and a lasso regression of the instrumented treatment variable, i.e., a regression of the predicted values of \code{d}, on controls \code{x}. 
#'
#' @param object an object of class \code{rlassoIV}, usually a result of a call \code{rlassoIV} with options \code{select.X=TRUE} and \code{select.Z=TRUE}.
#' @param selection.matrix if TRUE, a selection matrix is returned that indicates the selected variables from each first stage regression.
#' Default is set to FALSE. See section on details for more information.
#' @param complete general option of the function \code{coef}.
#' @param ... further arguments passed to function coef.
#' @return Coefficients obtained from \code{rlassoIV} by default. If option \code{selection.matrix} is \code{TRUE}, a list is returned with final coefficients, a matrix \code{selection.matrix}, and a matrix \code{selection.matrixZ}: 
#' \code{selection.matrix} contains the selection index for the lasso regression of \code{y} on \code{x} (first column) and the lasso regression of the predicted values of \code{d} on \code{x}
#' together with the union of these indizes.
#' \code{selection.matrixZ} contains the selection index from the first-stage lasso regression of \code{d} on \code{z} and \code{x}. 
#' @export
#' @rdname coef.rlassoIV
#' @examples 
#' \dontrun{
#' data(EminentDomain)
#' z <- EminentDomain$logGDP$z # instruments
#' x <- EminentDomain$logGDP$x # exogenous variables
#' y <- EminentDomain$logGDP$y # outcome varialbe
#' d <- EminentDomain$logGDP$d # treatment / endogenous variable
#' lasso.IV = rlassoIV(x=x, d=d, y=y, z=z, select.X=TRUE, select.Z=TRUE) 
#' coef(lasso.IV) # default behavior
#' coef(lasso.IV, selection.matrix = T) # print selection matrix
#' }
coef.rlassoIV <-  function(object, complete = TRUE, selection.matrix = FALSE, ...) {
  
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
    
      # selection w.r.t. 
      matZ <- object$selection.matrixZ
      dmatZ2 <- dim(matZ)[2]
      Zrnames <- rownames(matZ)
      matZ <- cbind(matZ, as.logical(apply(matZ, 1, sum)))
      colnames(matZ)[dim(matZ)[2]] <- "global.Z"
      matZ <- rbind(matZ, apply(matZ, 2, sum, na.rm = TRUE))
      matZ <- apply(matZ, 2, function(x) gsub(1, "x", x))
      matZ <- apply(matZ, 2, function(x) gsub(0, ".", x))
      # mat[is.na(matZ)] <- "-"
      rownames(matZ) <- c(Zrnames, "sum")
      
      if (complete) {
        coef <- list(cf = cf, selection.matrix = mat, 
                     selection.matrixZ = matZ)
        return(coef)
      }
      
      else {
        coef <- list(cf = cf[!is.na(cf)], selection.matrix = mat, 
                     selection.matrixZ = matZ)
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

############################################################################################

#' @rdname rlassoIV
#' @export

rlassoIVmult <- function(x, d, y, z, select.Z = TRUE, select.X = TRUE, 
                         ...) {
  # browser()
  d <- as.matrix(d)
  if (is.null(colnames(d))) 
    colnames(d) <- paste("d", 1:ncol(d), sep = "")
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:ncol(x), sep = "")
  if (is.null(colnames(z)) & !is.null(z)) 
    colnames(z) <- paste("z", 1:ncol(z), sep = "")
  
  if (select.Z == FALSE & select.X == FALSE) {
    res <- tsls(x=x, d=d, y=y, z=z, homoscedastic = FALSE, ...)
    return(res)
  }
  
  if (select.Z == TRUE & select.X == FALSE) {
    res <- rlassoIVselectZ(x, d, y, z, ...)
    return(res)
  }
  
  if (select.Z == FALSE & select.X == TRUE) {
    res <- rlassoIVselectX(x, d, y, z, ...)
    return(res)
  }
  
  if (select.Z == TRUE & select.X == TRUE) {
    d <- as.matrix(d)
    n <- dim(x)[1]
    d <- as.matrix(d)
    kd <- dim(d)[2]
    Z <- cbind(z, x)
    if (is.null(colnames(d))) 
      colnames(d) <- paste("d", 1:kd, sep = "")
    
    lasso.y.x <- rlasso(x, y, ...)
    Yr <- lasso.y.x$residuals
    Drhat <- NULL
    Zrhat <- NULL
    for (i in 1:kd) {
      lasso.d.x <- rlasso(d[, i] ~ x, ...)
      lasso.d.zx <- rlasso(d[, i] ~ Z, ...)
      if (sum(lasso.d.zx$index) == 0) {
        Drhat <- cbind(Drhat, d[, i] - mean(d[, i]))
        Zrhat <- cbind(Zrhat, d[, i] - mean(d[, i]))
        next
      }
      ind.dzx <- lasso.d.zx$index
      PZ <- Z[, ind.dzx, drop = FALSE] %*% MASS::ginv(t(Z[, ind.dzx, 
                                                          drop = FALSE]) %*% Z[, ind.dzx, drop = FALSE]) %*% t(Z[, 
                                                                                                                 ind.dzx, drop = FALSE]) %*% d[, i, drop = FALSE]
      lasso.PZ.x <- rlasso(PZ ~ x, ...)
      ind.PZx <- lasso.PZ.x$index
      Dr <- d[, i] - x[, ind.PZx, drop = FALSE] %*% MASS::ginv(t(x[, 
                                                                   ind.PZx, drop = FALSE]) %*% x[, ind.PZx, drop = FALSE]) %*% 
        t(x[, ind.PZx, drop = FALSE]) %*% PZ
      Zr <- lasso.PZ.x$residuals
      Drhat <- cbind(Drhat, Dr)
      Zrhat <- cbind(Zrhat, Zr)
    }
    result <- tsls(y = Yr, d = Drhat, x = NULL, z = Zrhat, homoscedastic = FALSE)
    coef <- as.vector(result$coefficient)
    se <- sqrt(diag(result$vcov))
    names(coef) <- names(se) <- colnames(d)
    res <- list(coefficients = coef, se = se, vcov = result$vcov, call = match.call(), 
                samplesize = n)
    class(res) <- "rlassoIV"
    return(res)
  }
}