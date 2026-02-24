source('R/01_estimate_dfb_evd.R')
source('R/05_rmaxdfbeta.R')

rdfbeta <- function(
  draw,
  n,
  nS = 1,
  xdist = \(n) rnorm(n),
  rdist = \(n) rnorm(n)
) {
  X <- matrix(xdist(draw * nS), ncol = nS)
  R <- matrix(rdist(draw * nS), ncol = nS)

  if (identical(xdist, \(n) rnorm(n))) {
    D <- rchisq(draw, n - nS)
  } else {
    D <- replicate(draw, sum(xdist(n - nS)^2))
  }

  dfb <- rowSums(R * X) / D
  return(dfb)
}

n <- c(25, 50, 100)
dfb <- sapply(n, \(i) rdfbeta(1e7, i, 1))
maxdfb <- sapply(n, \(i) replicate(1e5, rmaxdfbeta(i)))

pdf("paper/plots/fig1.pdf", width = 12, height = 5.5)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# --- Left plot: Delta_i ---
plot(
  density(dfb[, 3], n = 2000),
  xlim = c(0, 0.3),
  lwd = 2,
  main = expression(bold(Delta[i])),
  xlab = "",
  ylab = "",
  bty = "n"
)
grid(lty = 3, col = "lightgray")
lines(density(dfb[, 3], n = 2000), lwd = 2)
lines(density(dfb[, 2], n = 2000), lty = 4, col = "goldenrod", lwd = 2)
lines(density(dfb[, 1], n = 2000), lty = 2, col = "cyan3", lwd = 2)
legend(
  "topright",
  title = "Observations",
  legend = c("N = 100", "N = 50", "N = 25"),
  lty = c(1, 4, 2),
  col = c("black", "goldenrod", "cyan3"),
  lwd = 2,
  bty = "n"
)

# --- Right plot: max Delta_i ---
plot(
  density(maxdfb[, 3], n = 2000),
  xlim = c(0, 0.3),
  lwd = 2,
  main = expression(bold("max " * Delta[i])),
  xlab = "",
  ylab = "",
  bty = "n"
)
grid(lty = 3, col = "lightgray")
lines(density(maxdfb[, 3], n = 2000), lwd = 2)
lines(density(maxdfb[, 2], n = 2000), lty = 4, col = "goldenrod", lwd = 2)
lines(density(maxdfb[, 1], n = 2000), lty = 2, col = "cyan3", lwd = 2)
legend(
  "topright",
  title = "Observations",
  legend = c("N = 100", "N = 50", "N = 25"),
  lty = c(1, 4, 2),
  col = c("black", "goldenrod", "cyan3"),
  lwd = 2,
  bty = "n"
)

dev.off()
