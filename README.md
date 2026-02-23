# Testing Most Influential Sets

R implementation accompanying the paper [*"Testing Most Influential Sets"*](https://openreview.net/forum?id=zJ2AJz6xLO) by [Lucas D. Konrad](https://lucaskonrad.eu) and [Nikolas Kuschnig](https://www.kuschnig.eu/).

## Overview

Small influential data subsets can dramatically impact model conclusions — a few data points may overturn key findings. 
While recent work attempts to identify these most influential sets, there is no principled method to tell when influence of the maximally influential set is excessive rather than expected under natural sampling variation.

This package implements the statistical framework proposed in the paper for **principled hypothesis testing of most influential sets** in linear least-squares regression. 
We derive the extreme value distribution (EVD) of maximal influence and use it as a null distribution. 
Two distinct regimes emerge depending on set size relative to the sample: 
the heavy-tailed **Fréchet** distribution for constant-size sets with heavy-tailed data, and the well-behaved **Gumbel** distribution for growing sets or light-tailed data.

The implementation provides:

- An exact closed-form influence formula for observation sets (Proposition 1 of the paper), avoiding model re-fitting.
- EVD estimation via the block-maxima method with bias correction, using a greedy algorithm to approximate maximum influence within each block.
- Hypothesis testing for excessive influence: given a candidate set, estimate the null EVD and compute a p-value.
- The Frisch-Waugh-Lovell (FWL) theorem to partial out nuisance regressors, reducing the multivariate problem to a univariate one.
- Simulation utilities for Monte Carlo study of theoretical maximum influence distributions.

## Installation

```r
# Install external dependencies
install.packages(c("devtools", "evd"))

# Install from GitHub
devtools::install_github("konradld/testing_mis")
```

## Function Reference

### Core Functions

#### `estimate_dfb_evd(y, X1, Xother, S, block_count = 20)`

Main testing function implementing the three-step procedure from Section 3.2 of the paper. 
Uses FWL to partial out `Xother`, estimates tail coefficients of X and R to choose between Gumbel and Fréchet families (Theorem 1 / Corollary 1), 
then fits a GEV distribution to $\Delta$ block-maxima to assess whether the influence of set `S` is statistically excessive.

---


#### `dfbeta.numeric(Y, X, S, i_Xcol = 1)`

Computes the exact set influence $\Delta(S)$ for a subset `S` in a no-intercept simple regression of `Y` on a column of `X`. This implements the closed-form formula from Proposition 1 of the paper in the univariate (FWL-residualized) case.

---


#### `dfb_bmx(X, R, S, block_count)`

Constructs block-maxima of $\Delta(S)$ for EVD parameter estimation (Step 2 of the procedure). Divides the data (excluding `S`) into blocks, then uses an adaptive greedy algorithm within each block to approximate the most influential set of size `|S|`. This yields one maximum influence draw per block, suitable for GEV fitting.

---


#### `rmaxdfbeta(n = 100, xdist = rnorm, rdist = rnorm, nS = 1)`

Simulates a single draw from the distribution of $\Delta_{\max}$ over sets of size `nS` in a random sample of size `n`. Used for the convergence simulations in Section 4.1 of the paper (e.g., verifying Fréchet vs. Gumbel behavior under Normal and t(5) distributions).

---


### Helper Functions

#### `fwl(y, X1, X2)`

Implements the Frisch-Waugh-Lovell projection: residualizes `y` and `X1` on `X2` via QR decomposition. This reduces the multivariate regression to an equivalent univariate problem (see footnote 2 of the paper), enabling all subsequent influence computations in the simple regression setting.

---


#### `make_blocks(X, block_size)`

Reshapes a vector `X` into a matrix with `block_size` rows. Remainder elements beyond a full block are dropped.


## Dependency Graph

```
estimate_dfb_evd
├── fwl
├── dfbeta.numeric
├── make_blocks
├── dfb_bmx
│   └── make_blocks
└── evd::fgev

rmaxdfbeta
└── (standalone)

dfbeta.numeric
└── (standalone)
```

## Quick Example

The testing procedure follows three steps from Section 3.2 of the paper: (1) choose the EVD family based on tail behavior, (2) estimate EVD parameters via block maxima, and (3) perform the hypothesis test.

```r
library(evd)

set.seed(42)
n <- 500
X1 <- rnorm(n)
X2 <- rnorm(n)
y  <- 2 * X1 + X2 + rnorm(n)

# Suppose observations 1:5 are suspected to be influential
S <- 1:5

# Steps 1-2: estimate the EVD (selects Gumbel vs Fréchet, fits parameters)
result <- estimate_dfb_evd(y, X1, Xother = cbind(1, X2), S, block_count = 20)

# GEV parameters for the DFBETA block-maxima distribution
result$params

# Observed DFBETA of the suspicious set
result$set_dfb

# Step 3: p-value — is this influence excessive?
1 - evd::pgev(abs(result$set_dfb),
              loc   = result$params["loc"],
              scale = result$params["scale"],
              shape = result$params["shape"])
```

## Citation

```bibtex
@misc{konrad_testing_2025,
  title  = {Testing Most Influential Sets},
  author = {Konrad, Lucas D. and Kuschnig, Nikolas},
  year   = {2025},
  note   = {arXiv:2510.20372 [stat.ML]},
  url    = {https://arxiv.org/abs/2510.20372}
}
```

## License

[i dont know - pls check]
