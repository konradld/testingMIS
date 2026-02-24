#' Random Maximum Influence Draw
#'
#' Generates randomDFBETA values by simulating data.
#'
#' @param n_draw Number of draws
#' @param n Sample size (default: 100)
#' @param n_set Number of most influential observations to select (default: 1)
#' @param x_dist Function for generating X values (default: standard normal)
#' @param r_dist Function for generating residuals (default: standard normal)
#'
#' @return Numeric value of the set influence in a random sample
#'
#' @export
#' @importFrom stats rnorm rchisq
rdfbeta <- function(
  n_draw = 10L,
  n = 100L,
  n_set = 1L,
  x_dist = \(n) rnorm(n),
  r_dist = \(n) rnorm(n)
) {
  X <- matrix(x_dist(n_draw * n_set), ncol = n_set)
  R <- matrix(r_dist(n_draw * n_set), ncol = n_set)

  if (identical(x_dist, \(n) rnorm(n))) {
    D <- rchisq(n_draw, n - n_set)
  } else {
    D <- replicate(n_draw, sum(x_dist(n - n_set)^2))
  }

  dfb <- rowSums(R * X) / D
  return(dfb)
}
