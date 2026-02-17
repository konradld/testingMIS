# ===============================================================================
#
#               Functions for computing and bootstrapping DFBETA
#
# ===============================================================================


#' @description
#' Calculates the DFBETA value for a set of observations (S) using a computationally
#' efficient approach that avoids refitting the regression model multiple times.
#'
#' @param Y Numeric vector of response values
#' @param X Matrix or data frame of predictor variables
#' @param S Vector of indices for observations whose influence is being assessed
#' @param i_Xcol Column index of X to compute DFBETA for (default: 1)
#'
#' @return A named numeric value representing the DFBETA for observations in S
dfbeta.numeric <- function(
    Y,
    X,
    S,
    i_Xcol = 1) {
  if (!is.matrix(X)) { X <- as.matrix(X) }
  if (ncol(X) > 1) { X <- X[, i_Xcol] } # Extract X-column

  # Calculate sums
  sum_xy_all <- sum(X * Y)
  sum_x2_all <- sum(X^2)
  sum_xy_S <- sum(X[S] * Y[S])
  sum_x2_S <- sum(X[S]^2)

  # Calculate the result
  dfbeta <- -((sum_xy_all - sum_xy_S) / (sum_x2_all - sum_x2_S) - sum_xy_all / sum_x2_all)
  return("dfbeta" = dfbeta)
}


#' @description
#' Performs the non-parametric bootstrapping process to generate empirical distributions of DFBETA values
#'
#' @param NBoot Number of bootstrap iterations (default: 1000)
#' @param Y Response variable vector
#' @param X Matrix of predictor variables
#' @param n Total number of observations
#' @param nS Number of observations to sample in each bootstrap iteration
#' @param do_wild_bs Logical; whether to bootstrap x and r independently (default: FALSE)
#' @param do_normalize_x Logical; whether to normalize X variables (default: FALSE)
#' @param do_normalize_r Logical; whether to normalize residuals (default: FALSE)
#' @param do_normalize_dfb Logical; whether to normalize DFBETA values (default: FALSE)
#'
#' @return Numeric vector of length NBoot containing bootstrap DFBETA values
dfbeta_boot_nonparam <- function(
    NBoot = 1000,
    Y,
    X,
    n = length(Y),
    nS = 1,
    do_wild_bs = FALSE,
    do_normalize_x = FALSE,
    do_normalize_r = FALSE,
    do_normalize_dfb = FALSE) {

  if (do_normalize_x) { X <- scale(X) }
  if (do_normalize_r) {
    model <- lm(Y ~ X - 1)
    sigma_val <- sigma(model)
    Y <- Y / sigma_val
  }

  if (do_wild_bs) {
    if (do_normalize_r) {
      R <- model$residuals
    } else {
      R <- lm(Y ~ X - 1)$residuals
    }
    X_idx <- replicate(NBoot, sample(n, n, replace = TRUE), simplify = FALSE)
    R_idx <- replicate(NBoot, sample(n, nS, replace = TRUE), simplify = FALSE)
    boot_samples <- vapply(1:NBoot, \(i) {
      DFB <- 0
      for (j in 1:nS) {
        DFB_part <- X[X_idx[[i]][j], ] * R[R_idx[[i]][j]] / (sum(X[X_idx[[i]], ]^2) - sum(X[X_idx[[i]][1:j], ]^2))
        DFB <- DFB + DFB_part
      }
      return(DFB)
    }, numeric(1))
  } else {
    group_idx <- replicate(NBoot, sample(n, n, replace = TRUE), simplify = FALSE)
    boot_samples <- vapply(group_idx, \(i) dfbeta(Y = Y[i, drop = FALSE], X = X[i, , drop = FALSE], S = 1:nS), numeric(1))
  }

  if (do_normalize_dfb) {
    boot_samples <- scale(boot_samples)
  }

  return(boot_samples)
}


#' @description
#' Performs parametric bootstrap by fitting a Gaussian model to data,
#' then sampling from appropriate normal distributions for X and errors.
#'
#' @param NBoot Number of bootstrap iterations (default: 1000)
#' @param Y Response variable vector
#' @param X Matrix or vector of predictor variables
#' @param n Total number of observations
#' @param nS Number of observations in set S to analyze for influence
#' @param xdist Function for generating X values (default: standard normal)
#' @param rdist Function for generating residual values (default: standard normal)
#' @param do_normalize_x Logical; whether to normalize X variables (default: FALSE)
#' @param do_normalize_r Logical; whether to normalize residuals (default: FALSE)
#' @param do_normalize_dfb Logical; whether to normalize DFBETA values (default: FALSE)
#'
#' @return Numeric vector of length NBoot containing bootstrap DFBETA samples
dfbeta_boot_param <- function(
    NBoot = 1000,
    Y,
    X,
    n = length(Y),
    nS = 1,
    xdist = \(n) rnorm(n),
    rdist = \(n) rnorm(n),
    do_normalize_x = FALSE,
    do_normalize_r = FALSE,
    do_normalize_dfb = FALSE) {

  if (do_normalize_x) {
    X <- scale(X)
    mean_X <- 0
    sd_X <- 1
  } else {
    mean_X <- mean(X)
    sd_X <- sd(X)
  }
  if (do_normalize_r) {
    sd_R <- 1
  } else {
    model <- lm(Y ~ X - 1)
    residuals <- residuals(model)
    sd_R <- sd(residuals)
  }

  # Define distribution functions
  xdist_scaled <- \(n) mean_X + sd_X * xdist(n)
  rdist_scaled <- \(n) sd_R * rdist(n)

  # Call rdfbeta function
  boot_samples <- rdfbeta(
    draw = NBoot, n = n, nS = nS,
    xdist = xdist_scaled, rdist = rdist_scaled,
    do_normalize_x = do_normalize_x,
    do_normalize_r = do_normalize_r,
    do_normalize_dfb = do_normalize_dfb
  )

  return(boot_samples)
}


#' @description
#' Computes theoretical distribution samples for influence analysis in regression models.
#' This function is designed to be used within the parametric bootstrap framework.
#'
#' @param draw Number of samples
#' @param Y Response variable vector
#' @param X Matrix or vector of predictor variables
#' @param n Sample size (number of observations, default: length(Y))
#' @param nS Number of observations in set S to analyze for influence (default: 1)
#' @param thrtk_dist Function for generating theoretical distribution values
#' @param do_normalize_x Logical; whether to normalize X variables (default: FALSE)
#' @param do_normalize_r Logical; whether to normalize residuals (default: FALSE)
#' @param do_normalize_dfb Logical; whether to normalize DFBETA values (default: FALSE)
#'
#' @return Theoretical samples from rdfbeta function
dfbeta_theoretical <- function(
    draw = 1000,
    Y = NULL,
    X = NULL,
    n = length(Y),
    nS = 1,
    thrtk_dist,
    do_normalize_x = FALSE,
    do_normalize_r = FALSE,
    do_normalize_dfb = FALSE) {

  thrtk_samples <- rdfbeta(
    draw = draw,
    n = n,
    nS = nS,
    xdist = thrtk_dist,
    do_normalize_x = do_normalize_x,
    do_normalize_r = do_normalize_r,
    do_normalize_dfb = do_normalize_dfb
  )

  return(thrtk_samples)
}


#' @description
#' Performs bootstrap analysis to estimate the distribution of DFBETA values,
#' which measure the influence of observations on regression coefficients.
#'
#' @param Y Response variable vector
#' @param X Matrix of predictor variables
#' @param S Indices of observations to analyze for influence
#' @param nS Number of observations in set S
#' @param excl_set_bs Whether to exclude set S in bootstrap: "no" (default), "yes", or "both"
#' @param bootstrap_method Method to use: "nonparametric" (default) or "parametric"
#' @param wild_bs whether to bootstrap x and r independently: "no" (default), "yes", or "both"
#' @param NBoot Number of bootstrap iterations (default: 1000)
#' @param thrtk_dist Theoretical distribution (optional)
#' @param do_normalize_x Logical; whether to normalize X variables (default: FALSE)
#' @param do_normalize_r Logical; whether to normalize residuals (default: FALSE)
#' @param do_normalize_dfb Logical; whether to normalize DFBETA values (default: FALSE)
#' @param i_Xcol Column index of X to focus on (default: 1)
#'
#' @return A list containing bootstrap results
dfbeta_boot_new <- function(
    Y,
    X,
    S,
    nS,
    excl_set_bs = "no",
    bootstrap_method = "nonparametric",
    wild_bs = "yes",
    NBoot = 1000,
    thrtk_dist = NULL,
    do_normalize_x = FALSE,
    do_normalize_r = FALSE,
    do_normalize_dfb = FALSE,
    i_Xcol = 1) {

  if (!is.matrix(X)) { X <- as.matrix(X) }
  if (ncol(X) > 1) { X <- X[, i_Xcol] } # Extract X-column

  if (!excl_set_bs %in% c("no", "yes", "both")) {
    stop("excl_set_bs must be one of: 'no', 'yes', 'both'")
  }

  if (!bootstrap_method %in% c("nonparametric", "parametric", "both")) {
    stop("bootstrap_method must be one of: 'nonparametric', 'parametric', 'both'")
  }

  if (!wild_bs %in% c("no", "yes", "both")) {
    stop("bootstrap_method must be one of: 'no', 'yes', 'both'")
  }

  # Handle exclusion for main dataset
  original_X <- X
  original_Y <- Y

  if (excl_set_bs == "yes") {
    X <- X[-S, , drop = FALSE]
    Y <- Y[-S]
  }

  # Initialize variables
  n <- NROW(X)
  do_thrtk_dist <- !is.null(thrtk_dist)

  # Process based on bootstrap method
  result_list <- list()
  if (do_thrtk_dist) {
    result_list$theoretical <- dfbeta_theoretical(
      draw = NBoot,
      n = n,
      nS = nS,
      thrtk_dist = thrtk_dist,
      do_normalize_x = do_normalize_x,
      do_normalize_r = do_normalize_r,
      do_normalize_dfb = do_normalize_dfb
    )
  }
  if ((bootstrap_method == "nonparametric" || bootstrap_method == "both") && (wild_bs == "no" || wild_bs == "both")) {
    result_list$non_parametric_calm <- dfbeta_boot_nonparam(
      NBoot = NBoot,
      Y = Y,
      X = X,
      n = n,
      nS = nS,
      do_wild_bs = FALSE,
      do_normalize_x = do_normalize_x,
      do_normalize_r = do_normalize_r,
      do_normalize_dfb = do_normalize_dfb
    )
  }
  if ((bootstrap_method == "nonparametric" || bootstrap_method == "both") && (wild_bs == "yes" || wild_bs == "both")) {
    result_list$non_parametric_wild <- dfbeta_boot_nonparam(
      NBoot = NBoot,
      Y = Y,
      X = X,
      n = n,
      nS = nS,
      do_wild_bs = TRUE,
      do_normalize_x = do_normalize_x,
      do_normalize_r = do_normalize_r,
      do_normalize_dfb = do_normalize_dfb
    )
  }
  if (bootstrap_method == "parametric" || bootstrap_method == "both") {
    result_list$parametric <- dfbeta_boot_param(
      NBoot = NBoot,
      Y = Y,
      X = X,
      n = n,
      nS = nS,
      xdist = \(n) rnorm(n),
      rdist = \(n) rnorm(n),
      do_normalize_x = do_normalize_x,
      do_normalize_r = do_normalize_r,
      do_normalize_dfb = do_normalize_dfb
    )
  }

  # Run bootstrap under the alternative if requested
  if (excl_set_bs == "both") {
    X_excluded <- original_X[-S, , drop = FALSE]
    Y_excluded <- original_Y[-S, drop = FALSE]

    if (do_thrtk_dist) {
      result_list$theoretical_x <- dfbeta_theoretical(
        draw = NBoot,
        n = n - nS,
        nS = nS,
        thrtk_dist = thrtk_dist,
        do_normalize_x = do_normalize_x,
        do_normalize_r = do_normalize_r,
        do_normalize_dfb = do_normalize_dfb
      )
    }
    if ((bootstrap_method == "nonparametric" || bootstrap_method == "both") && (wild_bs == "no" || wild_bs == "both")) {
      result_list$non_parametric_calm_x <- dfbeta_boot_nonparam(
        NBoot = NBoot,
        Y = Y_excluded,
        X = X_excluded,
        n = n - nS,
        nS = nS,
        do_wild_bs = FALSE,
        do_normalize_x = do_normalize_x,
        do_normalize_r = do_normalize_r,
        do_normalize_dfb = do_normalize_dfb
      )
    }
    if ((bootstrap_method == "nonparametric" || bootstrap_method == "both") && (wild_bs == "yes" || wild_bs == "both")) {
      result_list$non_parametric_wild_x <- dfbeta_boot_nonparam(
        NBoot = NBoot,
        Y = Y_excluded,
        X = X_excluded,
        n = n - nS,
        nS = nS,
        do_wild_bs = TRUE,
        do_normalize_x = do_normalize_x,
        do_normalize_r = do_normalize_r,
        do_normalize_dfb = do_normalize_dfb
      )
    }
    if (bootstrap_method == "parametric" || bootstrap_method == "both") {
      result_list$parametric_x <- dfbeta_boot_param(
        NBoot = NBoot,
        Y = Y_excluded,
        X = X_excluded,
        n = n - nS,
        nS = nS,
        xdist = \(n) rnorm(n),
        rdist = \(n) rnorm(n),
        do_normalize_x = do_normalize_x,
        do_normalize_r = do_normalize_r,
        do_normalize_dfb = do_normalize_dfb
      )
    }
  }

  cat("Bootstrap analysis complete (method: ", bootstrap_method, ")\n", sep = "")
  return(result_list)
}


#' @description
#' Simulates a distribution of DFBETA or max(DFBETA) values using an optimized algorithm.
#'
#' @param draw Number of simulation iterations
#' @param n Sample size
#' @param nS Number of considered influential observatoins (default: 1)
#' @param xdist Function to generate X distribution (default: standard normal)
#' @param rdist Function to generate R (residuals) distribution (default: standard normal)
#' @param do_normalize_x Whether to normalize X variables (default: FALSE)
#' @param do_normalize_r Whether to normalize residuals (default: FALSE)
#' @param do_normalize_dfb Whether to normalize final DFBETA values (default: FALSE)
#'
#' @return A numeric vector of length `draw` containing (maximum) DFBETA values
#' from each simulation iteration.
rdfbeta <- function(
    draw, n, nS = 1,
    xdist = \(n) rnorm(n), rdist = \(n) rnorm(n),
    do_normalize_x = FALSE, do_normalize_r = FALSE, do_normalize_dfb = FALSE) {
  X <- matrix(xdist(draw * nS), ncol = nS)
  R <- matrix(rdist(draw * nS), ncol = nS)

  if (do_normalize_x) { X <- scale(X) }
  if (do_normalize_r) { R <- scale(R) }

  if (identical(xdist, \(n) rnorm(n))) {
    D <- rchisq(draw, n - nS)
  } else {
    D <- replicate(draw, sum(xdist(n - nS)^2))
  }

  dfb <- rowSums(R * X) / D
  if (do_normalize_dfb) { dfb <- scale(dfb) }

  return(dfb)
}

#===============================================================================
#               Functions for testing Maximum Infl. Set
#===============================================================================

#' @description
#' Estimates DFBETA parameters using extreme value distribution (EVD) fitting.
#' Uses the Frisch-Waugh-Lovell theorem to marginalize out other predictors 
#' and fits a generalized extreme value distribution to block maxima.
#' 
#' @param y Response variable vector
#' @param X1 Primary predictor variable of interest
#' @param Xother Matrix of other predictor variables to be marginalized out
#' @param S Vector of indices for observations whose influence is being assessed
#' @param block_count Number of blocks to divide the data into
#'
#' @return Numeric vector containing the estimated EVD parameters (location, scale, shape)
estimate_dfb_evd <- function(y, X1, Xother, S, block_count=20) {
  require(evd)
  
  fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
  
  Y <- fwl_vars[,1]
  X <- fwl_vars[,2]
  sumX2 <- sum(X^2)
  
  fwl_lm <- lm(Y~X-1)
  R <- fwl_lm$residuals
  
  Sdfb <- dfbeta(Y,X,S)
  
  # check if X and R are sufficiently small tailed
  bm_X <- apply(blocks(X[-S],length(X[-S])%/%block_count),2,max)
  bm_R <- apply(blocks(R[-S],length(R[-S])%/%block_count),2,max)
  
  x_evd <- evd::fgev(bm_X)
  r_evd <- evd::fgev(bm_R)
  
  x_ci <- c(x_evd$estimate['shape'] + 1.96*x_evd$std.err['shape'],
            x_evd$estimate['shape'] - 1.96*x_evd$std.err['shape'])
  
  r_ci <- c(r_evd$estimate['shape'] + 1.96*r_evd$std.err['shape'],
            r_evd$estimate['shape'] - 1.96*r_evd$std.err['shape'])
  
  # get block maxima
  Delta_bmx  <- abs(dfb_bmx(X, R, S = S, block_count))
  
  is_x_frechet <- all(x_ci > 0)
  is_r_frechet <- all(r_ci > 0)
  
  # if either is frechet, the larger significant one, 0 otherwise
  tail_coef <- max(is_x_frechet * x_evd$estimate['shape'],is_r_frechet * r_evd$estimate['shape']) 
  
  cat(
    sprintf(
      '\nX: Shape = %.4f \t| Pr(>x) = %.4f\nR: Shape = %.4f \t| Pr(>x) = %.4f\n',
      x_evd$estimate['shape'],
      pnorm(-x_evd$estimate['shape'], 0, x_evd$std.err['shape']),
      r_evd$estimate['shape'],
      pnorm(-r_evd$estimate['shape'], 0, r_evd$std.err['shape'])
    )
  )
  
  fit_evd_bm <- evd::fgev(Delta_bmx, shape = tail_coef)
  return(list('params'=fit_evd_bm$estimate, 'set_dfb' = Sdfb, 'block_maxima' = Delta_bmx))
}

#' @description
#' Detects DFBETA block-maxima for extreme value analysis. Divides data into blocks
#' and iteratively selects the most extreme DFBETA values within each block to
#' create a sample suitable for EVD fitting.
#' 
#' @param X Vector of predictor values
#' @param R Vector of residual values
#' @param S Vector of indices for the influential observation set
#' @param block_count Number of blocks to divide the data into
#'
#' @return Numeric vector of block maxima DFBETA values
dfb_bmx <- function(X, R, S, block_count) {
  # find out which tail we need
  Sdfb <- sum(X[S] * R[S]) / sum(X[-S]^2)
  if(Sdfb > 0) {
    which_tail <- 1
    which.extr <- \(x) which.max(x)
  } else if(Sdfb < 0) {
    which_tail <- -1
    which.extr <- \(x) which.min(x)
  } else {
    stop('dfbeta of S is exactly zero')
  }
  
  # delete infl. set from inference set
  X <- X[-S]
  R <- R[-S]
  
  # aux. variables and splitting into blocks
  sumX2 <- sum(X^2)
  nS  <- length(S)
  N0  <- length(X)
  block_size <- N0%/%block_count
  Xbl <- blocks(X, block_size) # wasteful, I dont allow for different sized blocks, so leftovers are dropped
  Rbl <- blocks(R, block_size) # wasteful, I dont allow for different sized blocks, so leftovers are dropped
  
  # Initialize list to store selection matrices
  dfb_list <- list()
  m_list <- list()
  
  # First selection
  dfb_list[[1]] <- sapply(1:block_count, \(i) 
                          Xbl[,i] * Rbl[,i] / (sumX2 - Xbl[,i]^2))
  m_list[[1]] <- apply(dfb_list[[1]], 2, which.extr)
  
  # Subsequent selections
  if (nS > 1) {
    for (s in 2:nS) {
      prev_selected <- m_list[1:(s-1)]
      
      dfb_list[[s]] <- sapply(1:block_count, \(i) {
        # Get indices of previously selected items
        selected_idx <- sapply(prev_selected, \(m) m[i])
        # Calculate for remaining items
        remaining <- setdiff(1:block_size, selected_idx)
        
        dfb_vals <- Xbl[remaining,i] * Rbl[remaining,i] / 
          (sumX2 - sum(Xbl[selected_idx,i]^2) - Xbl[remaining,i]^2)
        
        # Return with original indices
        result <- numeric(block_size)
        result[remaining] <- dfb_vals
        result[selected_idx] <- - which_tail * Inf  # Ensure already selected aren't picked again
        result
      })
      
      m_list[[s]] <- apply(dfb_list[[s]], 2, which.extr)
    }
  }
  
  # Combine all selected indices and compute final result
  all_selected <- do.call(rbind, m_list)
  
  sapply(1:block_count, \(i) 
         sum(Xbl[all_selected[,i], i] * Rbl[all_selected[,i], i]) / 
           (sumX2 - sum(Xbl[all_selected[,i], i]^2)))
}

#' @description
#' Generates random maximum DFBETA values by simulating data and iteratively
#' selecting the most influential observations. Used for theoretical distribution
#' analysis of maximum influence measures.
#' 
#' @param n Sample size (default: 100)
#' @param xdist Function for generating X values (default: standard normal)
#' @param rdist Function for generating residual values (default: standard normal)
#' @param nS Number of most influential observations to select (default: 1)
#'
#' @return Numeric value representing the maximum influence for a set in a random sample
rmaxdfbeta <- function(n = 100, xdist = \(n)rnorm(n), rdist = \(n)rnorm(n), nS = 1) {
  # Generate random data
  x <- xdist(n)
  r <- rdist(n)
  
  # x <- xdist(n%/%block_count)
  # r <- rdist(n%/%block_count)
  
  x_sq <- sum(x^2)
  # x_sq <- sum(xdist(n)^2)
  
  # Initialize vectors to store selected indices
  selected_indices <- integer(nS)
  
  # Keep track of remaining indices
  remaining <- seq_len(n)
  
  # Cumulative sum of squared x values for selected indices
  x_sq_selected <- 0
  
  # Iteratively find the nS indices that maximize the influence measure
  for (s in seq_len(nS)) {
    # Calculate the influence measure for remaining indices
    # The denominator adjusts for previously selected x values
    influence <- (x[remaining] * r[remaining]) / 
      (x_sq - x_sq_selected - x[remaining]^2)
    
    # Find which remaining index maximizes the influence
    max_idx_in_remaining <- which.max(influence)
    
    # Store the actual index in the original data
    selected_indices[s] <- remaining[max_idx_in_remaining]
    
    # Update the cumulative sum of squared x values for selected indices
    x_sq_selected <- x_sq_selected + x[selected_indices[s]]^2
    
    # Remove the selected index from remaining indices
    remaining <- remaining[-max_idx_in_remaining]
  }
  
  # Calculate and return the final measure
  # Sum of (x * r) for selected indices divided by adjusted x_sq
  numerator <- sum((x * r)[selected_indices])
  denominator <- x_sq - x_sq_selected
  
  return(numerator / denominator)
}

rmaxdfbeta_brute <- function(n = 100, xdist = \(n)rnorm(n), rdist = \(n)rnorm(n), nS = 1) {
  # Generate random data
  x <- xdist(n)
  r <- rdist(n)
  sumX2 <- sum(x^2)
  
  max(combn(1:n, nS, \(i) sum(x[i] * r[i]) / (sumX2 - sum(x[i]^2))))
}

# Frisch-Waugh-Lovell Theorem to marginalize out X2
fwl <- function(y, X1, X2) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  XtX2_inv <- solve(crossprod(X2))
  Yp <- y - X2 %*% (XtX2_inv %*% crossprod(X2, y))
  Xp <- X1 - X2 %*% (XtX2_inv %*% crossprod(X2, X1))
  return(cbind('Yp' = Yp, 'Xp' = Xp))
}

# Helper function for cutting data into equal sized blocks (wastefully)
blocks <- function(X, block_size) {
  n_complete_blocks <- floor(length(X) / block_size)
  
  start_indices <- seq(1, by = block_size, length.out = n_complete_blocks)
  sapply(start_indices, \(i) X[i:(i + block_size - 1)])
}