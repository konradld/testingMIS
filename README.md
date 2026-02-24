# Testing Most Influential Sets

R package accompanying the paper [*"Testing Most Influential Sets"*](https://openreview.net/forum?id=zJ2AJz6xLO) by [Lucas D. Konrad](https://lucaskonrad.eu) and [Nikolas Kuschnig](https://www.kuschnig.eu/).

## Overview

Small subsets of observations can dramatically change regression estimates -- a handful of data points may overturn key findings. While methods exist to identify these most influential sets, there has been no principled way to determine whether the observed influence is excessive or simply expected under normal sampling variation.

This package implements a statistical framework for **hypothesis testing of most influential sets** in linear least-squares regression. The approach derives the extreme value distribution (EVD) of maximal influence and uses it as a null distribution. Two regimes emerge depending on set size and tail behavior:

- **Frechet** -- for constant-size sets with heavy-tailed data.
- **Gumbel** -- for growing sets or light-tailed data.

The implementation provides:

- An exact closed-form influence formula for observation sets (Proposition 1), avoiding model re-fitting.
- EVD estimation via the block-maxima method with a greedy algorithm to approximate maximum influence within each block.
- Hypothesis testing: given a candidate set, estimate the null EVD and compute a p-value.
- The Frisch-Waugh-Lovell (FWL) theorem to partial out nuisance regressors, reducing the problem to a univariate one.
- Simulation utilities for Monte Carlo study of maximum influence distributions.

## Installation

```r
# Install from GitHub
devtools::install_github("konradld/testing_mis")
```

## Quick Example

The testing procedure follows three steps (Section 3.2 of the paper): (1) classify the EVD family based on tail behavior, (2) estimate EVD parameters via block maxima, and (3) compute a p-value for the candidate set's influence.

```r
library("testingMIS")

set.seed(42)
n <- 500L
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- 2 * x1 + x2 + rnorm(n)

# Suppose observations 1:5 are suspected to be influential
S <- 1:5

# Steps 1-2: estimate the EVD (classifies Gumbel vs Frechet, fits parameters)
result <- estimate_dfb_evd(y, x = x1, Z = cbind(1, x2), set = S, block_count = 20)

# GEV parameters for the DFBETA block-maxima distribution
result$params

# Observed DFBETA of the suspicious set
result$set_dfb

# Step 3: p-value -- is this influence excessive?
1 - evd::pgev(
  abs(result$set_dfb),
  loc = result$params["loc"],
  scale = result$params["scale"],
  shape = result$params["shape"]
)
```

## Function Reference

### `estimate_dfb_evd(y, x, Z, set, block_count = 20, verbose = TRUE)`

Main testing function implementing the three-step procedure from Section 3.2. Uses FWL to partial out `Z`, estimates tail coefficients of X and R to classify Gumbel vs Frechet (Theorem 1 / Corollary 1), then fits a GEV distribution to block maxima. Returns a list with:

- `$params` -- named vector of GEV parameters (`loc`, `scale`, `shape`).
- `$set_dfb` -- observed DFBETA of the candidate set.
- `$block_maxima` -- the block-maxima values used for fitting.

Set `verbose = FALSE` to suppress diagnostic output.

---

### `dfbeta_numeric(y, X, set, col_X = 1)`

Computes the exact set influence for a subset `set` in a no-intercept regression of `y` on column `col_X` of `X`. Implements the closed-form formula from Proposition 1 in the univariate (FWL-residualized) case.

An S3 method `dfbeta.numeric` is also registered, allowing dispatch via `stats::dfbeta()` on numeric vectors.

---

### `dfb_bmx(X, R, S, block_count)`

Constructs block maxima of set influence for EVD parameter estimation (Step 2). Divides the data (excluding `S`) into blocks, then uses a greedy algorithm within each block to approximate the most influential set of size `|S|`. Returns one maximum influence value per block, suitable for GEV fitting.

---

### `rmaxdfbeta(n = 100L, n_set = 1L, x_dist = rnorm, r_dist = rnorm)`

Simulates a single draw from the distribution of the maximum DFBETA over sets of size `n_set` in a random sample of size `n`. The `x_dist` and `r_dist` arguments are functions taking a single argument (sample size) that generate predictor and residual values, respectively. Used for convergence simulations (e.g., verifying Frechet vs Gumbel behavior).

---

### `rdfbeta(n_draw = 10L, n = 100L, n_set = 1L, x_dist = rnorm, r_dist = rnorm)`

Generates `n_draw` random DFBETA values by simulation. Each draw computes the set influence for a random set of size `n_set` in a sample of size `n`. Useful for Monte Carlo study of the influence distribution (as opposed to the *maximum* influence distribution targeted by `rmaxdfbeta`).

---
## Replication

The `paper/` directory contains scripts to replicate all figures, tables, and applications from the paper.

**Simulations:**

- `10_simulate_fig1.R` -- Figure 1 (EVD convergence)
- `11_simulate_fig2.R` -- Figure 2
- `12_simulate_fig3_tabA1.R` -- Figure 3 and Table A1
- `13_simulate_figA1.R` -- Figure A1

**Applications:**

- `20_application_rugged.R` -- Rugged terrain and economic development
- `21_application_sparrows.R` -- Sparrows dataset
- `22_application_law.R` -- Law School Admissions
- `23_application_adult.R` -- UCI Adult income
- `24_application_housing.R` -- Housing prices
- `25_application_communities.R` -- Communities and Crime

Output is stored in `paper/plots/` (PDF figures) and `paper/tables/` (LaTeX tables).

### Datasets

The `data/` directory contains datasets used in the paper's empirical applications:

| File | Description | Application script |
|------|-------------|--------------------|
| `adult.data` | UCI Adult income dataset | `paper/23_application_adult.R` |
| `communities.data` | UCI Communities and Crime dataset | `paper/25_application_communities.R` |
| `lsa.RData` | Law School Admissions data | `paper/22_application_law.R` |
| `SparrowsElphick.txt` | Elphick's Sparrows dataset | `paper/21_application_sparrows.R` |

## Citation

```bibtex
@inproceedings{konrad_kuschnig_2026,
  title     = {Testing Most Influential Sets},
  author    = {Konrad, Lucas D. and Kuschnig, Nikolas},
  booktitle = {Proceedings of the International Conference on Learning Representations (ICLR)},
  year      = {2026},
  note      = {arXiv:2510.20372 [stat.ML]},
  url       = {https://arxiv.org/abs/2510.20372}
}
```
