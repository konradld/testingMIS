# Linear Regression Model with Influential Set Analysis
# Publication-ready version

rm(list = ls())
try(dev.off())

set.seed(1234)

# Load required libraries
library(evd)

# Source custom functions (update path as needed)
source('R/4_bootstrap-dfb.R')

# Helper function for blocking
blocks <- function(X, block_size) {
  sapply(seq(1, length(X) - block_size + 1, block_size), 
         \(i) X[i:(i + block_size - 1)])
}

# ============================================
# Data Generation
# ============================================

# Generate base data
n <- 100
x <- seq(1, 10, length.out = n)
mx <- mean(x)
x <- x - mx

# True relationship: y = 3*x + error
y <- 3*x + rnorm(n, 0, 2)

# Add influential point(s)
influential_x <- c(11) - mx  # High leverage point
influential_y <- c(23.5)      # Outlier value

# Combine data
x_all <- c(influential_x, x)
y_all <- c(influential_y, y)

nS <- length(influential_x)
S <- 1:nS

# Create data frame
data <- data.frame(
  x = x_all,
  y = y_all,
  type = c(rep("Influential", nS), rep("Regular", n))
)

# ============================================
# Model Fitting
# ============================================

# Fit models
model_full <- lm(y ~ x - 1, data = data)
data_clean <- data[data$type == "Regular", ]
model_clean <- lm(y ~ x - 1, data = data_clean)

# Calculate DFBETA for influential set
Sdfb <- dfbeta(y_all, x_all, S = S)

# ============================================
# Block Maximum Analysis
# ============================================

X <- x
sumX2 <- sum(X^2)
R <- model_full$residuals
N <- n
block_count <- 10

# Calculate block maxima
D_bsmx <- {
  Xbl <- blocks(X, N/block_count)
  Rbl <- blocks(R, N/block_count)
  
  dfb1 <- sapply(1:block_count, 
                 \(i) Xbl[,i] * Rbl[,i] / (sumX2 - Xbl[,i]^2))
  m1 <- apply(dfb1, 2, which.max)
  
  dfb2 <- sapply(1:block_count, 
                 \(i) Xbl[-m1[i],i] * Rbl[-m1[i],i] / 
                   (sumX2 - Xbl[m1[i],i]^2 - Xbl[-m1[i],i]^2))
  m2 <- apply(dfb2, 2, which.max)
  
  m2 <- m2 + (m1 < m2)  # Fix index
  sapply(1:block_count, 
         \(i) sum(Xbl[c(m1[i], m2[i]),i] * Rbl[c(m1[i], m2[i]),i]) / 
           (sumX2 - sum(Xbl[c(m1[i], m2[i]),i]^2)))
}

# Fit Gumbel distribution
fit_evd_bs <- evd::fgev(D_bsmx, shape = 0)
param_bs <- fit_evd_bs$estimate

# Calculate quantiles
ex90qant <- evd::qgumbel(0.90, 
                         loc = param_bs[1] + param_bs[2] * log(block_count), 
                         scale = param_bs[2])
ex95qant <- evd::qgumbel(0.95, 
                         loc = param_bs[1] + param_bs[2] * log(block_count), 
                         scale = param_bs[2])
ex99qant <- evd::qgumbel(0.99, 
                         loc = param_bs[1] + param_bs[2] * log(block_count), 
                         scale = param_bs[2])

# Calculate p-value
pval <- 1 - evd::pgumbel(Sdfb, 
                         loc = param_bs[1] + param_bs[2] * log(block_count), 
                         scale = param_bs[2])

# ============================================
# Publication-Ready Plotting
# ============================================

# Set up publication-quality plot parameters
par(mfrow = c(1, 2), 
    mar = c(4.5, 4.5, 3, 1), 
    mgp = c(2.5, 0.7, 0),
    las = 1,
    cex.lab = 1.1,
    cex.axis = 0.9,
    cex.main = 1.2)

# Color palette (colorblind-friendly)
col_influential <- "#D55E00"  # Orange-red
col_regular <- "#0072B2"      # Blue
col_full_fit <- "#CC79A7"     # Reddish purple
col_clean_fit <- "#009E73"    # Bluish green
col_confidence <- c("#999999", "#666666", "#333333")  # Grays for CI bands

# ---- Panel A: Regression with Influential Points ----
plot(x_all, y_all, 
     pch = c(17, rep(16, n)),  # Triangle for influential, circle for regular
     cex = c(1.5, rep(1, n)),  # Larger size for influential point
     col = c(col_influential, rep(col_regular, n)),
     ylim = c(-50, 50),
     xlim = range(x_all) * 1.05,
     main = expression(bold("(a) Effect of Influential Observation")),
     xlab = expression(italic(x)),
     ylab = expression(italic(y)),
     panel.first = grid(col = "gray90", lty = 1))

# Add regression lines
abline(model_clean, col = col_clean_fit, lwd = 2.5)
abline(model_full, col = col_full_fit, lwd = 2.5, lty = 2)

# Add confidence bands for single influential point
if(nS == 1) {
  par_coef <- model_clean$coefficients
  sumX2clean <- sum((data_clean$x)^2)
  
  # 90% confidence band
  curve(par_coef * x + sumX2clean * ex90qant / x, 
        add = TRUE, col = col_confidence[1], lwd = 1.5, lty = 3)
  curve(par_coef * x - sumX2clean * ex90qant / x, 
        add = TRUE, col = col_confidence[1], lwd = 1.5, lty = 3)
  
  # 95% confidence band
  curve(par_coef * x + sumX2clean * ex95qant / x, 
        add = TRUE, col = col_confidence[2], lwd = 1.5, lty = 3)
  curve(par_coef * x - sumX2clean * ex95qant / x, 
        add = TRUE, col = col_confidence[2], lwd = 1.5, lty = 3)
  
  # 99% confidence band
  curve(par_coef * x + sumX2clean * ex99qant / x, 
        add = TRUE, col = col_confidence[3], lwd = 1.5, lty = 3)
  curve(par_coef * x - sumX2clean * ex99qant / x, 
        add = TRUE, col = col_confidence[3], lwd = 1.5, lty = 3)
}

# Add legend
legend("topleft", 
       legend = c("Regular observations", 
                  "Influential observation",
                  "Fit without influential", 
                  "Fit with influential",
                  "90%, 95%, 99% bands"),
       pch = c(16, 17, NA, NA, NA),
       col = c(col_regular, col_influential, 
               col_clean_fit, col_full_fit, col_confidence[2]),
       lty = c(NA, NA, 1, 2, 3),
       lwd = c(NA, NA, 2.5, 2.5, 1.5),
       pt.cex = c(1, 1.5, NA, NA, NA),
       bg = "white",
       box.col = "gray80",
       cex = 0.85)

# ---- Panel B: Gumbel Distribution Fit ----
# Calculate densities
dens_empirical <- density(D_bsmx, bw = "SJ")
x_range <- seq(0, max(Sdfb * 1.4, max(D_bsmx)), length.out = 200)
y_gumbel <- evd::dgumbel(x_range, 
                         loc = param_bs[1], 
                         scale = param_bs[2])

y_gumbel_cor <- evd::dgumbel(x_range, 
                             loc = param_bs[1] + param_bs[2] * log(block_count), 
                             scale = param_bs[2])

# Create empty plot with proper limits
plot(NULL, NULL,
     xlim = c(0, max(Sdfb, max(D_bsmx))*1.4),
     ylim = c(0, max(max(dens_empirical$y), max(y_gumbel)) * 1.1),
     main = expression(bold("(b) Extreme Value Analysis")),
     xlab = "DFBETA statistic",
     ylab = "Density",
     panel.first = grid(col = "gray90", lty = 1))

# Add histogram in background
hist(D_bsmx, 
     breaks = "FD",
     freq = FALSE,
     col = rgb(0, 114, 178, alpha = 30, maxColorValue = 255),
     border = NA,
     add = TRUE)

# Add Gumbel fit
lines(x_range, y_gumbel, 
      col = col_clean_fit, 
      lwd = 2.5, 
      lty = 2)

# Add Gumbel fit
lines(x_range, y_gumbel_cor, 
      col = "black", 
      lwd = 2.5, 
      lty = 1)

# Add observed statistic
abline(v = Sdfb, 
       col = col_influential, 
       lwd = 2.5, 
       lty = 1)

# Add shaded area for p-value
x_shade <- seq(Sdfb, max(x_range), length.out = 100)
y_shade <- evd::dgumbel(x_shade, 
                        loc = param_bs[1] + param_bs[2] * log(block_count), 
                        scale = param_bs[2])
polygon(c(x_shade, rev(x_shade)), 
        c(y_shade, rep(0, length(y_shade))),
        col = rgb(213, 94, 0, alpha = 50, maxColorValue = 255),
        border = NA)

# Add text annotation for p-value
text(Sdfb * 1.05, max(y_gumbel) * 0.9,
     paste0("p = ", formatC(pval, format = "f", digits = 3)),
     pos = 4,
     cex = 1.1,
     font = 2)

# Add legend
legend("topright",
       legend = c("Block maxima",
                  "Corrected Gumbel",
                  "Uncorrected Gumbel fit",
                  "Observed statistic"),
       fill = c(rgb(0, 114, 178, alpha = 30, maxColorValue = 255),
                NA, NA, NA),
       border = c(NA, NA, NA, NA),
       lty = c(NA, 1, 2, 1),
       col = c(NA, "black", col_clean_fit, col_influential),
       lwd = c(NA, 2.5, 2.5, 2.5),
       bg = "white",
       box.col = "gray80",
       cex = 0.85)

# ============================================
# Optional: Save to file
# ============================================

# Uncomment to save as PDF
dev.copy(pdf, "paper/plots/fig2.pdf", width = 12, height = 5.5)
dev.off()

# Print summary statistics
cat("\n=== Analysis Summary ===\n")
cat("Number of observations:", n, "\n")
cat("Number of influential points:", nS, "\n")
cat("Regression coefficient (clean):", round(model_clean$coefficients, 3), "\n")
cat("Regression coefficient (full):", round(model_full$coefficients, 3), "\n")
cat("DFBETA statistic:", round(Sdfb, 3), "\n")
cat("P-value:", formatC(pval, format = "f", digits = 4), "\n")
cat("90% quantile:", round(ex90qant, 3), "\n")
cat("95% quantile:", round(ex95qant, 3), "\n")
cat("99% quantile:", round(ex99qant, 3), "\n")