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
  source('R/04_make_blocks.R')
  
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