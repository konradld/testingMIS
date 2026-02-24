#' Set DFBETA for numeric inputs
#'
#' Calculates the DFBETA for a set of observations (S)
#'
#' @param Y Numeric vector of response values
#' @param X Matrix or data frame of predictor variables
#' @param S Vector of indices for observations whose influence is being assessed
#' @param i_Xcol Column index of X to compute DFBETA for (default: 1)
#'
#' @return A numeric value representing the DFBETA for observations in S
#'
#' @export
dfbeta_numeric <- function(Y, X, S, i_Xcol = 1) {
  X <- as.numeric(as.matrix(X)[, i_Xcol])
  sum_xy_all <- sum(X * Y)
  sum_x2_all <- sum(X^2)
  c(
    dfbeta = -((sum_xy_all - sum(X[S] * Y[S])) /
      (sum_x2_all - sum(X[S]^2)) -
      sum_xy_all / sum_x2_all)
  )
}

#' @export
#' @method dfbeta numeric
dfbeta.numeric <- function(model, ...) {
  dfbeta_numeric(Y = model, ...)
}

#' @importFrom stats dfbeta
NULL
