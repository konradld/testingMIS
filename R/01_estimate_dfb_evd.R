#' Estimate EVD of Set Influence
#'
#' Estimate the extreme value distribution (EVD) parameters of set influence.
#' Uses the Frisch-Waugh-Lovell theorem to marginalize out other predictors
#' and fits a generalized extreme value distribution to block maxima.
#'
#' @param y Response variable vector
#' @param X1 Primary predictor variable of interest
#' @param Xother Matrix of other predictor variables to be marginalized out
#' @param S Vector of indices for observations whose influence is being assessed
#' @param block_count Number of blocks to divide the data into
#' @param verbose Logical indicating whether to print results
#'
#' @return List containing the estimated EVD parameters, set DFBETA, and block maxima
#'
#' @export
#' @importFrom stats lm residuals pnorm rnorm
#' @importFrom evd fgev
#'
estimate_dfb_evd <- function(
  y,
  X1,
  Xother,
  S,
  block_count = 20,
  verbose = TRUE
) {
  fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)

  Y <- fwl_vars[, 1]
  X <- fwl_vars[, 2]

  fwl_lm <- lm(Y ~ X - 1)
  R <- fwl_lm |> residuals()

  Sdfb <- dfbeta_numeric(Y, X, S)

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

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\nX: Shape = %.4f \t| Pr(>x) = %.4f\nR: Shape = %.4f \t| Pr(>x) = %.4f\n",
      x_evd$estimate["shape"],
      pnorm(-x_evd$estimate["shape"], 0, x_evd$std.err["shape"]),
      r_evd$estimate["shape"],
      pnorm(-r_evd$estimate["shape"], 0, r_evd$std.err["shape"])
    ))
  }

  # Get block maxima and fit EVD

  Delta_bmx <- abs(dfb_bmx(X, R, S = S, block_count))
  fit_evd_bm <- evd::fgev(Delta_bmx, shape = tail_coef)
  fit_evd_bm$estimate["shape"] <- tail_coef

  list(params = fit_evd_bm$estimate, set_dfb = Sdfb, block_maxima = Delta_bmx)
}
