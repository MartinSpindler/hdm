#' Backend for rlassoATEAutoDML
#'
#' Data backend for auto debiased ATE. Implements the preprocessing steps required for estimation of the ATE using auto debiased ATE.
#' 
#' 
#' @param x \code{character} specifying exogenous variables
#' @param d \code{character} specifying treatment variable (binary)
#' @param y \code{character} specifying outcome variable / dependent variable
#' @param data \code{data.frame} or \code{data.table} containing the raw input data.
#' @param dict a dictionary (\code{function})
#' function of (D,X) with D being a numeric and X a matrix, which returns a numeric vector. Default is \code{dict <- function(d, x) return(c(1, d, x))}. The dictionary will be applied on an individual level, i.e., being executed for each i in 1, ..., n with i representing a row in the \code{data}.
#' @return an object of class \code{DataAutoDebiasedATE} with the following entries
#' \itemize{
#'   \item{\code{data_prep} - A copy of the original \code{data} (containing only variables as provided through \code{x}, \code{d} and \code{y})}
#'   \item{\code{Y} - Y as passed as input} 
#'   \item{\code{X} - Processed X (dictionary applied to input data)}
#'   \item{\code{p} - Number of columns of preprocessed matrix \code{X}}
#'   \item{\code{M} - Matrix M as of ... TODO }
#'   \item{\code{X_prelim} - Version of matrix X used for a preliminary estimator}
#'   \item{\code{M_prelim} - Version of matrix M used for a preliminary estimator}
#'   \item{\code{p_prelim_out} - Number of columns of preprocessed matrix \code{X_prelim}}
#' }
#' 
#' 
#' @export
DataATEAutoDML <- function(x, d, y, data, dict = NULL) {
  # TODO: Alternatively provide input as matrices or provide wrapper
  checkmate::check_character(x)
  checkmate::check_character(d, max.len = 1)
  checkmate::check_character(y, max.len = 1)
  
  if (!is.null(dict)) {
    checkmate::check_class(dict, "function")
  } else {
    dict <- default_dict_ATE
  }
  if (all(class(data) == "data.frame")) {
    data = data.table::data.table(data)
  }
  checkmate::assert_class(data, "data.table")
  n <- nrow(data)
  all_vars <- c(x, d, y)
  data_prep <- data.table::copy(data[, ..all_vars])
  
  X_raw <- as.matrix(data_prep[, ..x])
  D <- data_prep[[d]]
  Y <- data_prep[[y]]
  
  # Compute X: Result from dictionary applied to D and Xs
  X1 <- dict(D[1], X_raw[1,])
  p <- length(X1)
  checkmate::check_numeric(X1)
  X <- matrix(nrow = n, ncol = p)
  for (i in seq_len(n)) {
    X[i, ] <- dict(D[i], X_raw[i, ])
  }
  # Precompute objects for preliminary estimator (based on random choice)
  if (p > 60) {
    p_prelim <- min(ceiling(p/40), p)
  } else {
    p_prelim <- min(ceiling(p/4), p)
  }
  x_vars_prelim <- x[1:p_prelim]
  X_raw_prelim <- X_raw[, x_vars_prelim, drop = FALSE]
  p_prelim_out <- length(dict(D[1], X_raw_prelim[1, ]))
  X_prelim <- matrix(nrow = n, ncol = p_prelim_out)

  for (i in seq_len(n)) {
   X_prelim[i, ] <- dict(D[i], X_raw_prelim[i, ]) 
  }
  # Compute M: Difference for dictionary applied for D = 1 and D = 0 for all obs
  M <- matrix(nrow = n, ncol = p)
  
  for (i in seq_len(n)) {
    M[i, ] <- dict(1, X_raw[i, ]) - dict(0, X_raw[i, ])
  }
  # Precompute M for preliminary estimator
  M_prelim <- matrix(nrow = n, ncol = p_prelim_out)
  for (i in seq_len(n)) {
    M_prelim[i, ] <- dict(1, X_raw_prelim[i, ]) - dict(0, X_raw_prelim[i, ])
  }
  
  obj = list("data_prep" = data_prep, "M" = M,
             "Y" = Y, "X" = X, "p" = p, "n" = n,
             "M_prelim" = M_prelim, "X_prelim" = X_prelim,
             "p_prelim_out" = p_prelim_out)
  class(obj) <- "DataATEAutoDML"
  return(obj)
}


#' Backend for rlassoAPDAutoDML
#'
#' Data backend for auto-debiased estimation of average partial derivatives 
#' (continuous treatment variable). Implements the preprocessing steps required 
#' for auto-debiased estimation of the APD. If interested in using panel data, 
#' user must also enter "unit" and "time" variables as well. If using cross section
#' these fields can be left empty (\code{NULL}). 
#' 
#' 
#' @param x \code{character} specifying exogenous variables
#' @param d \code{character} specifying treatment variable (continuous)
#' @param y \code{character} specifying outcome variable / dependent variable.
#' @param x_manual \code{character} specifying variables that should manually be added to the model. These variables are excluded from the the internal preprocessing w.r.t. to polynomials and interactions. The corresponding data will be adjusted during derivative calculations. Variables should not include interactions with variable \code{d}.
#' @param data \code{data.frame} or \code{data.table} with data set
#' @param poly_order Polynomial order.
#' @param interactions \code{logical} specifying whether two-way interactions of \code{x} and \code{d} should be constructed. Default is \code{true}.
#' @param unit \code{character} specifying unit variable of data in panel context, if not panel data omit (default \code{NULL}).
#' @param time \code{character} specifying time variable of data in panel context, if not panel data omit (default \code{NULL}).
#' @return an object of class \code{DataAPDAutoDML} with the following entries
#' \itemize{
#'   \item{\code{data_prep} - \code{data.table} object containing the preprocessed data}
#'   \item{\code{Y} - Y as passed as input} 
#'   \item{\code{X} - Processed X (dictionary applied to input data)}
#'   \item{\code{p} - Number of columns of preprocessed matrix \code{X}}
#'   \item{\code{M} - Matrix M as of ... TODO }
#' }
#'   
#' @export
DataAPDAutoDML <- function(x, d, y, x_manual = NULL, data, poly_order = 3,
                           interactions = TRUE, unit = NULL, 
                           time = NULL) {
  checkmate::check_character(x)
  checkmate::check_character(d, max.len = 1)
  checkmate::check_character(y, max.len = 1)
  
  if (!is.null(x_manual)) {
    checkmate::check_character(x_manual)
  }
  
  if (all(class(data) == "data.frame")) {
    data = data.table::data.table(data)
  }
  checkmate::assert_class(data, "data.table")
  checkmate::assert_character(names(data), unique = TRUE)
  checkmate::check_numeric(poly_order, lower = 2, len = 1)
  
  # Check if columns in data have positive variance and drop otherwise
  const_cols <- names(which(apply(data, 2, function(x) var(x) < 10^(-5))))
  
  if (length(const_cols) > 0) {
    data[, c(const_cols):=NULL]
    x <- x[x != const_cols ]
    x_manual <- x_manual[x_manual != const_cols]
  }
  
  if(!is.null(unit) & !is.null(time)){
    checkmate::check_character(unit, max.len = 1)
    checkmate::check_character(time, max.len = 1)
    # panel data clean by sorting 
    setorderv(data, c(unit, time))
  }
  
  p <- length(x)
  n <- nrow(data)
  all_vars <- c(x, d, y, x_manual)
  
  data_prep <- data.table::copy(data[, ..all_vars])
  d_cols <- character(length = poly_order)
  d_cols[1] <- d
  
  # X columns for polynomials (excluding binary x variables)
  x_bin_indx <- apply(data_prep[, ..x], 2, function(x) { all(x %in% 0:1) })
  x_bin <- names(data_prep[, ..x])[x_bin_indx]
  x_notbin <- setdiff(x, x_bin)
  p_notbin <- length(x_notbin)
  
  if (p_notbin > 0) {
    x_cols <- character(length = (p_notbin*poly_order))
    x_cols[1:p_notbin] <- x_notbin
    dx_names <- c(d, x_notbin)
  } else {
    x_cols <- x
    dx_names <- d
  }
  
  # construct polynomials
  for (pol_ord in 2:poly_order) {
    d_cols_to_poly <- paste0(d, "_", pol_ord)
    d_cols[pol_ord] <- d_cols_to_poly
    new_cols_to_poly <- c(d_cols_to_poly)
    
    if (p_notbin > 0) {
      x_cols_to_poly <- paste0(x_notbin, "_", pol_ord)
      x_cols[(((pol_ord - 1)*p_notbin + 1) : ((pol_ord - 1)*p_notbin + p_notbin))] <- x_cols_to_poly
      new_cols_to_poly <- c(d_cols_to_poly, x_cols_to_poly)
    }
    data_prep[, paste0(new_cols_to_poly) := lapply(.SD, function(x) { x^pol_ord }), .SDcols = dx_names] 
  }
  # generate interactions between covars x
  if (interactions) {
    cols_int_x_x <- paste0("int_", combn(x,2, FUN=paste, collapse='_')) 
    comb_covars <- combn(x, 2)
    for (col in 1:length(comb_covars[1,])){
      col_name_new <- paste0("int_", comb_covars[1, col], "_", 
                             comb_covars[2, col] )
      data_prep[[col_name_new]] <- data_prep[[comb_covars[1, col]]]*data_prep[[comb_covars[2, col]]]
      if(length(unique(data_prep[[col_name_new]])) <= 1){
        data_prep[[col_name_new]] <- NULL
      }
    }
    cols_int_x_x <- intersect(colnames(data_prep), cols_int_x_x) # vector of final interactions 
    x_cols <- unique(c(x_cols, x_bin, cols_int_x_x))
  } else {
    x_cols <- unique(c(x_cols, x_bin))
  }
  
  # generate interactions for d with x
  if (interactions) {
    for (this_x in x) {
      interaction_col_name <- paste0("int_", d, "_", this_x)
      data_prep[[interaction_col_name]] <- data_prep[[d]] * data_prep[[this_x]]
    }
  }
  
  # construct matrices M: manual derivative
  cols_with_covars <- setdiff(names(data_prep), y)
  
  M <- data.table::copy(data_prep[, ..cols_with_covars ])
  M[[d]] <- rep(1, n)
  
  for (pol_ord in 2:poly_order) {
    this_d <- d_cols[pol_ord]
    prev_d <- d_cols[pol_ord - 1]
    M[[this_d]] <-  pol_ord * data_prep[[prev_d]]
  }
  
  # make all all columns that are not a function of D have deriv 0
  M[, (x_cols) := lapply(.SD, function(x) { rep(0, n) }), .SDcols = x_cols]
  
  if (!is.null(x_manual)) {
    M[, (x_manual) := lapply(.SD, function(x) { rep(0, n) }), .SDcols = x_manual]
  }
  
  if (interactions) {
    # derivative of X interacted with D is just X 
    for (this_x in x){
      this_int <- paste0("int_", d, "_", this_x)
      M[[this_int]] <- data_prep[[this_x]]
    }

  }
    
  # Normalize M
  # save the SD of the X matrix - we use these to normalize the M 
  basis_sd <- GMCM:::colSds(as.matrix(data_prep[,..cols_with_covars])) 
  M <- t(t(M)/basis_sd)
  
  # normalize all covars in data_prep
  data_prep[, (cols_with_covars) := lapply(.SD, function(x) { scale(x) }), .SDcols = cols_with_covars]

  # if panel data take first difference 
  if(!is.null(unit) & !is.null(time)){
    #handle data_prep first 
    ## add back unit and time 
    data_prep$unit <- data[[unit]]
    data_prep$time <- data[[time]]
    all_vars <- c(cols_with_covars, y)
    data_prep <- data_prep[, (all_vars  ) := lapply(.SD, function(v) {v - shift(v)}), .(unit), .SDcols = all_vars ]
    first_time_period <- min(data_prep$time)
    idx_first_time_period <- data_prep[time == first_time_period, which = TRUE]
    data_prep <- data_prep[- idx_first_time_period, ] # removing first time period bc of first diff
    # now also remove first time period from M, but we dont first diff M 
    M <- M[-idx_first_time_period,]
  }
  
  # export Y and X as vector & matrix
  Y <- data_prep[[y]]
  X <- as.matrix(data_prep[, ..cols_with_covars])
  M <- M
  p <- ncol(X)
  n <- nrow(X) 
  # Output
  # data_prep: Preprocessed data table
  # Y, X, M: Matrices
  obj = list("data_prep" = data_prep, "M" = M,
             "Y" = Y, "X" = X, "p" = p, "n" = n)
  class(obj) <- "DataAPDAutoDML"
  return(obj)
}
