#' @description
#' Generates random maximum DFBETA values by simulating data and iteratively selecting
#' the most influential observations in a greedy manner for theoretical distribution
#' analysis of maximum influence measures.
#'
#' @param n Sample size (default: 100)
#' @param xdist Function for generating X values (default: standard normal)
#' @param rdist Function for generating residual values (default: standard normal)
#' @param nS Number of most influential observations to select (default: 1)
#'
#' @return Numeric value representing the maximum influence for a set in a random sample
rmaxdfbeta <- function(n = 100, xdist = rnorm, rdist = rnorm, nS = 1) {
  x <- xdist(n)
  r <- rdist(n)
  x_sq <- sum(x^2)
  
  selected <- integer(nS)
  remaining <- seq_len(n)
  x_sq_selected <- 0
  
  for (s in seq_len(nS)) {
    influence <- (x[remaining] * r[remaining]) /
      (x_sq - x_sq_selected - x[remaining]^2)
    
    best <- which.max(influence)
    selected[s] <- remaining[best]
    x_sq_selected <- x_sq_selected + x[selected[s]]^2
    remaining <- remaining[-best]
  }
  
  sum((x * r)[selected]) / (x_sq - x_sq_selected)
}