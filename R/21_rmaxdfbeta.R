#' Random Maximum Influence Draw
#'
#' Generates random maximum DFBETA values by simulating data and iteratively
#' selecting the most influential observations in a greedy manner for
#' theoretical distribution analysis of maximum influence measures.
#'
#' @param n Sample size (default: 100)
#' @param n_set Number of most influential observations to select (default: 1)
#' @param x_dist Function for generating X values (default: standard normal)
#' @param r_dist Function for generating residuals (default: standard normal)
#'
#' @return Numeric value of the maximum set influence in a random sample
#'
#' @export
rmaxdfbeta <- function(
  n = 100L,
  n_set = 1L,
  x_dist = \(n) rnorm(n),
  r_dist = \(n) rnorm(n)
) {
  x <- x_dist(n)
  r <- r_dist(n)
  x_sq <- sum(x^2)

  selected <- integer(n_set)
  remaining <- seq_len(n)
  x_sq_selected <- 0

  for (s in seq_len(n_set)) {
    influence <- (x[remaining] * r[remaining]) /
      (x_sq - x_sq_selected - x[remaining]^2)

    best <- which.max(influence)
    selected[s] <- remaining[best]
    x_sq_selected <- x_sq_selected + x[selected[s]]^2
    remaining <- remaining[-best]
  }

  sum((x * r)[selected]) / (x_sq - x_sq_selected)
}
