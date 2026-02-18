rm(list = ls())
try(dev.off())
set.seed(1234)

library(evd)
source('R/4_bootstrap-dfb.R')

# ── Data ─────────────────────────────────────────────────────────────────────
n <- 100
x <- scale(seq(1, 10, length.out = n), scale = FALSE)[,1]
y <- 3*x + rnorm(n, 0, 2)

S  <- 1L
x_all <- c(11 - mean(seq(1, 10, length.out = n)), x)
y_all <- c(23.81, y)
type  <- c("Influential", rep("Regular", n))

data       <- data.frame(x = x_all, y = y_all, type)
model_full <- lm(y ~ x - 1, data = data)
model_clean <- lm(y ~ x - 1, data = data, subset = type == "Regular")
Sdfb       <- dfbeta.numeric(y_all, x_all, S = S)

# ── Block Maxima ──────────────────────────────────────────────────────────────
block_count <- 10L
R_all  <- model_full$residuals
R      <- R_all[-1]

sumX2  <- sum(x^2)
Xbl    <- matrix(x, nrow = n/block_count)
Rbl    <- matrix(R, nrow = n/block_count)

D_bsmx <- sapply(1:block_count, \(i) {
  d1 <- Xbl[,i] * Rbl[,i] / (sumX2 - Xbl[,i]^2)
  m1 <- which.max(d1)
  d2 <- Xbl[-m1,i] * Rbl[-m1,i] / (sumX2 - Xbl[m1,i]^2 - Xbl[-m1,i]^2)
  m2 <- which.max(d2) + (which.max(d1) < which.max(d2))
  sum(Xbl[c(m1,m2),i] * Rbl[c(m1,m2),i]) / (sumX2 - sum(Xbl[c(m1,m2),i]^2))
})

# ── EVD Fit ───────────────────────────────────────────────────────────────────
fit    <- fgev(D_bsmx, shape = 0)$estimate
loc_adj <- fit[1] + fit[2] * log(block_count)
quants <- qgumbel(c(.90, .95, .99), loc = loc_adj, scale = fit[2])
pval   <- 1 - pgumbel(Sdfb, loc = loc_adj, scale = fit[2])

# ── Plot ──────────────────────────────────────────────────────────────────────
pal <- list(inf = "#D55E00", reg = "#0072B2", full = "#CC79A7",
            clean = "#009E73", ci = c("#999999","#666666","#333333"))

par(mfrow = c(1,2), mar = c(4.5,4.5,3,1), mgp = c(2.5,0.7,0),
    las = 1, cex.lab = 1.1, cex.axis = 0.9, cex.main = 1.2)

# Panel A
plot(x_all, y_all,
     pch = c(17, rep(16, n)), cex = c(1.5, rep(1, n)),
     col = c(pal$inf, rep(pal$reg, n)),
     ylim = c(-50, 50), xlim = range(x_all) * 1.05,
     main = expression(bold("(a) Effect of Influential Observation")),
     xlab = expression(italic(x)), ylab = expression(italic(y)),
     panel.first = grid(col = "gray90", lty = 1))

abline(model_clean, col = pal$clean, lwd = 2.5)
abline(model_full,  col = pal$full,  lwd = 2.5, lty = 2)

b0        <- model_clean$coefficients
sumX2cl   <- sum(model_clean$model$x^2)
for (i in seq_along(quants)) {
  curve(b0*x + sumX2cl*quants[i]/x, add=TRUE, col=pal$ci[i], lwd=1.5, lty=3)
  curve(b0*x - sumX2cl*quants[i]/x, add=TRUE, col=pal$ci[i], lwd=1.5, lty=3)
}

legend("bottomright", bg = "white", box.col = "gray80", cex = 0.85,
       legend = c("Regular","Influential","Fit without","Fit with","90/95/99% bands"),
       pch = c(16,17,NA,NA,NA), pt.cex = c(1,1.5,NA,NA,NA),
       col = c(pal$reg, pal$inf, pal$clean, pal$full, pal$ci[2]),
       lty = c(NA,NA,1,2,3), lwd = c(NA,NA,2.5,2.5,1.5))

# Panel B
x_seq <- seq(0, max(Sdfb, max(D_bsmx)) * 1.4, length.out = 200)
yG    <- dgumbel(x_seq, loc_adj, fit[2])
dens  <- density(D_bsmx, bw = "SJ")

plot(NULL, xlim = range(x_seq), ylim = c(0, max(dens$y, yG) * 1.1),
     main = expression(bold("(b) Extreme Value Analysis")),
     xlab = "DFBETA statistic", ylab = "Density",
     panel.first = grid(col = "gray90", lty = 1))

hist(D_bsmx, breaks = "FD", freq = FALSE, add = TRUE,
     col = rgb(0,114,178,30,maxColorValue=255), border = NA)
lines(x_seq, dgumbel(x_seq, fit[1], fit[2]), col = pal$clean, lwd = 2.5, lty = 2)
lines(x_seq, yG, lwd = 2.5)
abline(v = Sdfb, col = pal$inf, lwd = 2.5)

x_sh <- seq(Sdfb, max(x_seq), length.out = 100)
polygon(c(x_sh, rev(x_sh)), c(dgumbel(x_sh, loc_adj, fit[2]), rep(0, 100)),
        col = rgb(213,94,0,50,maxColorValue=255), border = NA)
text(Sdfb*1.05, max(yG)*0.9, sprintf("p = %.3f", pval),
     pos=4, cex=1.1, font=2)

legend("topright", bg="white", box.col="gray80", cex=0.85,
       legend = c("Block maxima","Corrected Gumbel","Uncorrected Gumbel","Observed"),
       fill = c(rgb(0,114,178,30,maxColorValue=255),NA,NA,NA),
       border = NA, lty = c(NA,1,2,1),
       col = c(NA,"black",pal$clean,pal$inf), lwd = c(NA,2.5,2.5,2.5))

# ── Output ────────────────────────────────────────────────────────────────────
dev.copy(pdf, "paper/plots/fig2.pdf", width = 12, height = 5.5); dev.off()

cat(sprintf("\n=== Analysis Summary ===
n = %d | Influential pts = 1
β (clean) = %.3f | β (full) = %.3f
DFBETA = %.3f | p-value = %.4f
Quantiles: 90%% = %.3f | 95%% = %.3f | 99%% = %.3f\n",
            n, model_clean$coefficients, model_full$coefficients,
            Sdfb, pval, quants[1], quants[2], quants[3]))