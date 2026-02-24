#' Set DFBETA for numeric inputs
#'
#' Calculates the DFBETA for a set of observations (S)
#'
#' @param y Numeric vector of response values
#' @param X Matrix or data frame of predictor variables
#' @param set Vector of indices for the set whose influence is being assessed
#' @param col_X Column index of X to compute DFBETA for (default: 1)
#'
#' @return A numeric value representing the DFBETA for observations in S
#'
#' @export
dfbeta_numeric <- function(y, X, set, col_X = 1) {
  X <- as.matrix(X)[, col_X]
  sum_xy_all <- sum(X * y)
  sum_xsq_all <- sum(X^2)
  c(
    dfbeta = -((sum_xy_all - sum(X[set] * y[set])) /
      (sum_xsq_all - sum(X[set]^2)) -
      sum_xy_all / sum_xsq_all)
  )
}

#' @export
#' @method dfbeta numeric
dfbeta.numeric <- function(model, ...) {
  dfbeta_numeric(y = model, ...)
}

#' @importFrom stats dfbeta
NULL
