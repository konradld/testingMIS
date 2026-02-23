# ===============================================================================
#
#               Functions for computing and bootstrapping DFBETA
#
# ===============================================================================


#' @description
#' Calculates the DFBETA value for a set of observations (S)
#'
#' @param Y Numeric vector of response values
#' @param X Matrix or data frame of predictor variables
#' @param S Vector of indices for observations whose influence is being assessed
#' @param i_Xcol Column index of X to compute DFBETA for (default: 1)
#'
#' @return A numeric value representing the DFBETA for observations in S
dfbeta.numeric <- function(Y, X, S, i_Xcol = 1) {
  X <- as.numeric(as.matrix(X)[, i_Xcol])
  sum_xy_all <- sum(X * Y)
  sum_x2_all <- sum(X^2)
  c(dfbeta = -((sum_xy_all - sum(X[S] * Y[S])) / (sum_x2_all - sum(X[S]^2)) - sum_xy_all / sum_x2_all))
}


# ===============================================================================
#               Functions for testing Maximum Infl. Set
# ===============================================================================

#' @description
#' Estimates extreme value distribution (EVD) of DFBETA
#' Uses the Frisch-Waugh-Lovell theorem to marginalize out other predictors
#' and fits a generalized extreme value distribution to block maxima.
#'
#' @param y Response variable vector
#' @param X1 Primary predictor variable of interest
#' @param Xother Matrix of other predictor variables to be marginalized out
#' @param S Vector of indices for observations whose influence is being assessed
#' @param block_count Number of blocks to divide the data into
#'
#' @return List containing the estimated EVD parameters, set DFBETA, and block maxima
estimate_dfb_evd <- function(y, X1, Xother, S, block_count = 20) {
  fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
  
  Y <- fwl_vars[, 1]
  X <- fwl_vars[, 2]
  
  fwl_lm <- lm(Y ~ X - 1)
  R <- fwl_lm$residuals
  
  Sdfb <- dfbeta.numeric(Y, X, S)
  
  # Check if X and R are sufficiently small-tailed
  bm_X <- apply(make_blocks(X[-S], length(X[-S]) %/% block_count), 2, max)
  bm_R <- apply(make_blocks(R[-S], length(R[-S]) %/% block_count), 2, max)
  
  x_evd <- evd::fgev(bm_X)
  r_evd <- evd::fgev(bm_R)
  
  is_x_frechet <- x_evd$estimate["shape"] - 1.96 * x_evd$std.err["shape"] > 0
  is_r_frechet <- r_evd$estimate["shape"] - 1.96 * r_evd$std.err["shape"] > 0
  
  # If either is Frechet, take the larger significant shape; 0 otherwise
  tail_coef <- max(
    is_x_frechet * x_evd$estimate["shape"],
    is_r_frechet * r_evd$estimate["shape"]
  )
  
  cat(sprintf(
    "\nX: Shape = %.4f \t| Pr(>x) = %.4f\nR: Shape = %.4f \t| Pr(>x) = %.4f\n",
    x_evd$estimate["shape"],
    pnorm(-x_evd$estimate["shape"], 0, x_evd$std.err["shape"]),
    r_evd$estimate["shape"],
    pnorm(-r_evd$estimate["shape"], 0, r_evd$std.err["shape"])
  ))
  
  # Get block maxima and fit EVD
  
  Delta_bmx <- abs(dfb_bmx(X, R, S = S, block_count))
  fit_evd_bm <- evd::fgev(Delta_bmx, shape = tail_coef)
  
  list(params = fit_evd_bm$estimate, set_dfb = Sdfb, block_maxima = Delta_bmx)
}

#' @description
#' Detects DFBETA block-maxima for extreme value analysis. Divides data into blocks
#' and iteratively selects the most extreme DFBETA values greedily within each block to
#' create a sample suitable for EVD fitting.
#'
#' @param X Vector of predictor values
#' @param R Vector of residual values
#' @param S Vector of indices for the influential observation set
#' @param block_count Number of blocks to divide the data into
#'
#' @return Numeric vector of block maxima DFBETA values
dfb_bmx <- function(X, R, S, block_count) {
  sgn <- sign(sum(X[S] * R[S]))
  if (sgn == 0) stop("dfbeta of S is exactly zero")
  which.extr <- if (sgn > 0) which.max else which.min
  
  # Remove influential set from inference set
  X <- X[-S]
  R <- R[-S]
  
  sumX2 <- sum(X^2)
  nS <- length(S)
  block_size <- length(X) %/% block_count
  Xbl <- make_blocks(X, block_size)
  Rbl <- make_blocks(R, block_size)
  
  # First selection: compute per-element DFBETA within each block
  dfb_mat <- Xbl * Rbl / (sumX2 - Xbl^2)
  m_list <- list(apply(dfb_mat, 2, which.extr))
  
  # Subsequent selections
  if (nS > 1) {
    for (s in 2:nS) {
      dfb_mat <- vapply(seq_len(block_count), function(i) {
        selected_idx <- vapply(m_list, function(m) m[i], integer(1))
        remaining <- setdiff(seq_len(block_size), selected_idx)
        
        result <- rep(-sgn * Inf, block_size)
        result[remaining] <- Xbl[remaining, i] * Rbl[remaining, i] /
          (sumX2 - sum(Xbl[selected_idx, i]^2) - Xbl[remaining, i]^2)
        result
      }, numeric(block_size))
      
      m_list[[s]] <- apply(dfb_mat, 2, which.extr)
    }
  }
  
  # Combine all selected indices and compute final result
  all_selected <- do.call(rbind, m_list)
  
  vapply(seq_len(block_count), function(i) {
    idx <- all_selected[, i]
    sum(Xbl[idx, i] * Rbl[idx, i]) / (sumX2 - sum(Xbl[idx, i]^2))
  }, numeric(1))
}

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


# ===============================================================================
#               Helper functions
# ===============================================================================

#' Frisch-Waugh-Lovell projection: residualizes y and X1 on X2
fwl <- function(y, X1, X2) {
  Q <- qr(as.matrix(X2))
  cbind(Yp = y - qr.fitted(Q, y), Xp = as.matrix(X1) - qr.fitted(Q, as.matrix(X1)))
}

#' Splits a vector into equal-sized blocks (drops remainder)
make_blocks <- function(X, block_size) {
  n <- length(X) %/% block_size
  matrix(X[seq_len(n * block_size)], nrow = block_size, ncol = n)
}