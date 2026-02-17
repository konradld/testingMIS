rm(list = ls())
try(dev.off())
# set.seed(1234)
source('./00_code/dfbeta_funcs.R')

# Load required library
library(evd)

# Parameters
nS = 1
N = 5000 + nS
sx = sqrt(2)
sr = sqrt(2)
xdist <- function(n) rnorm(n, 0, sx)
rdist <- function(n) rnorm(n, 0, sr)

# xdist <- function(n) rt(n, 2)
# rdist <- function(n) rt(n, 2)
block_count = 30

# Number of densities to average
M <- 1000 # Number of replications for averaging

# Function to generate samples and fit both Gumbel methods
fit_both_gumbel_methods <- function() {
  success <- FALSE
  while (!success) {
    res <- tryCatch({
      maxdraws <- replicate(1000, rmaxdfbeta(n = N-nS, xdist = xdist, rdist = rdist, nS = nS))
      fit_evd <- fgev(maxdraws, std.err = FALSE, warn.inf = FALSE, shape = 0) # == change here for Frechet
      param <- fit_evd$estimate
      
      # Method 2: Block maxima
      X <- xdist(N)
      U <- rdist(N)
      Y <- 2 * X + U
      
      lm_mod <- lm(Y ~ X - 1)
      R <- resid(lm_mod)
      
      S <- 1:nS
      
      D_bmx <- abs(dfb_bmx(X, R, S, block_count)) # get BM
      fit_evd_bm <- fgev(D_bmx, std.err = FALSE, warn.inf = FALSE, shape = 0) # ====== change here for Frechet
      param_bm <- fit_evd_bm$estimate
      
      success <- TRUE  # Set flag to TRUE if no error occurs
      
      return(list(
        regular = param,
        block_maxima = param_bm
      ))
    }, error = function(e) {
      cat("Error occurred:", e$message, "\nRetrying...\n")
      NULL  # Return NULL on error
    })
  }
  return(res)
}

# Generate M sets of Gumbel parameters for both methods
cat("Generating", M, "samples and fitting Gumbel distributions (both methods)...\n")
results <- replicate(M, fit_both_gumbel_methods(), simplify = FALSE)

# Extract parameters for both methods
regular_params <- t(sapply(results, function(x) x$regular))
bm_params <- t(sapply(results, function(x) x$block_maxima))
ratios <- t(sapply(results, function(x) x$ratio))

# Calculate average parameters
mean_loc_regular <- mean(regular_params[, 1])
mean_scale_regular <- mean(regular_params[, 2])
mean_shape_regular <- ifelse(ncol(regular_params) == 3, mean(regular_params[, 3]), NA)
mean_loc_bm <- mean(bm_params[, 1])
mean_scale_bm <- mean(bm_params[, 2])
mean_loc_bm_cor <- mean(bm_params[, 1] + bm_params[, 2] * log(block_count))
mean_shape_bm <- ifelse(ncol(bm_params) == 3, mean(bm_params[, 3]), NA) 

cat("\n=== Regular Method ===\n")
cat("Average location:", mean_loc_regular, "\n")
cat("Average scale:", mean_scale_regular, "\n")

cat("\n=== Block Maxima Method ===\n")
cat("Average location:", mean_loc_bm, "\n")
cat("Average scale:", mean_scale_bm, "\n")

cat("\n=== Average Ratios ===\n")
cat("Location ratio (regular/bm):", mean(ratios[, 1]), "\n")
cat("Scale ratio (regular/bm):", mean(ratios[, 2]), "\n")

# Create sequence for plotting
x_seq <- seq(0, 0.01, length.out = 1000)

# Calculate average densities for both methods
avg_density_regular <- numeric(length(x_seq))
avg_density_bm <- numeric(length(x_seq))
avg_density_bm_adjusted <- numeric(length(x_seq))

for (i in 1:M) {
  # Regular method
  avg_density_regular <- avg_density_regular + 
    dgumbel(x_seq, loc = regular_params[i, 1], scale = regular_params[i, 2])
  
  # Block maxima - unadjusted
  avg_density_bm <- avg_density_bm + 
    dgumbel(x_seq, loc = bm_params[i, 1], scale = bm_params[i, 2])
  
  # Block maxima - adjusted with log(block_count)
  avg_density_bm_adjusted <- avg_density_bm_adjusted + 
    dgumbel(x_seq, 
            loc = bm_params[i, 1] + bm_params[i, 2] * log(block_count), 
            scale = bm_params[i, 2])
}

avg_density_regular <- avg_density_regular / M
avg_density_bm <- avg_density_bm / M
avg_density_bm_adjusted <- avg_density_bm_adjusted / M



# ============================================
# Plotting for Gumbel Distribution Comparison
# ============================================

# Set up publication-quality plot parameters

par(mfrow = c(1, 3 + (ncol(bm_params) == 3)), 
    mar = c(4.5, 4.5, 3, 2), 
    mgp = c(2.8, 0.8, 0),
    las = 1,
    cex.lab = 1.15,
    cex.axis = 0.95,
    cex.main = 1.2,
    font.lab = 1,
    font.main = 1)

# Professional color palette (colorblind-friendly)
col_regular <- "#0072B2"      # Blue
col_bm <- "#009E73"           # Orange-red  
col_bm_adjusted <- "#D55E00"  # Green
col_avg_param <- "#999999"    # Gray
col_hist <- "#A6D8F4"         # Light blue
col_hist2 <- "#F0E442"        # Yellow

# ---- Panel A: Density Comparison ----
plot(x_seq, avg_density_regular, 
     type = "l", 
     lwd = 2, 
     col = col_regular,
     main = expression(bold("(a) Average Gumbel Densities")),
     xlab = expression(Delta^max),
     ylab = "Density",
     ylim = c(0, max(c(avg_density_regular, avg_density_bm, avg_density_bm_adjusted)) * 1.05),
     # xlim = range(x_seq),
     xlim = c(0,0.003),
     ann = TRUE,       # turn off annotations
     axes = FALSE)

# axis(1, at = pretty(range(x_seq)))  
axis(1, at = c(0,0.001,0.002,0.003))  
axis(2, las = 1) 
# Add density lines
lines(x_seq, avg_density_bm, col = col_bm, lwd = 2, lty = 2)
lines(x_seq, avg_density_bm_adjusted, col = col_bm_adjusted, lwd = 2, lty = 4)

# Add densities using average parameters (thinner, dashed)
# lines(x_seq, dgumbel(x_seq, loc = mean_loc_regular, scale = mean_scale_regular),
#       col = col_regular, lwd = 1.5, lty = 2)
# lines(x_seq, dgumbel(x_seq, loc = mean_loc_bm, scale = mean_scale_bm),
#       col = col_bm, lwd = 1.5, lty = 2)
# lines(x_seq, dgumbel(x_seq, loc = mean_loc_bm + mean_scale_bm * log(block_count), 
#                      scale = mean_scale_bm),
#       col = col_bm_adjusted, lwd = 1.5, lty = 2)

# Professional legend
legend("topright", 
       legend = c(
         "  True",
         "  BM", 
         "  BM (adj.)"
       ),
       col = c(col_regular, col_bm, col_bm_adjusted),
       lwd = c(2, 2, 2),
       lty = c(1, 2, 4),
       bty = "n",
       cex = 1.0,
       y.intersp = 0.9,
       inset = 0.1,
       seg.len = 4)  # Add this line - default is usually 2

# Add subtle box around plot area
# box(lwd = 1.2, col = "gray40")

# ---- Panel B: Location Parameter Distribution ----
# Calculate adjusted locations
adjusted_locations <- bm_params[, 1] + bm_params[, 2] * log(block_count)

# Create histogram
h1 <- hist(adjusted_locations, 
           breaks = 25,
           plot = FALSE)

plot(h1, 
     main = expression(bold("(b) Location Parameter"~tilde(a))),
     xlab = expression(tilde(a)),
     ylab = "Frequency",
     col = col_hist,
     border = "white",
     lwd = 0.5,
     xlim = range(c(adjusted_locations, mean_loc_regular)) * c(0.95, 1.05)
)

# Add density overlay
d1 <- density(adjusted_locations, bw = "SJ")
par(new = TRUE)
plot(d1$x, d1$y * diff(h1$mids[1:2]) * length(adjusted_locations),
     type = "l",
     lwd = 2,
     col = "gray30",
     xlim = range(c(adjusted_locations, mean_loc_regular)) * c(0.95, 1.05),
     ylim = c(0, max(h1$counts)),
     axes = FALSE,
     xlab = "",
     ylab = "")

# Add vertical lines for means
abline(v = mean(adjusted_locations), col = col_bm_adjusted, lwd = 2, lty = 1)
abline(v = mean_loc_regular, col = col_regular, lwd = 2, lty = 2)

# Add mean values as text
text(mean(adjusted_locations), max(h1$counts) * 0.95, 
     paste0("BM: ", round(mean(adjusted_locations), 3)),
     pos = 2, cex = 0.9, font = 2, col = col_bm_adjusted)
text(mean_loc_regular, max(h1$counts) * 0.85,
     paste0("True: ", round(mean_loc_regular, 3)),
     pos = 4, cex = 0.9, font = 2, col = col_regular)

# box(lwd = 1.2, col = "gray40")

# ---- Panel C: Scale Parameter Distribution ----
h2 <- hist(bm_params[, 2], 
           breaks = 25,
           plot = FALSE)

plot(h2,
     main = expression(bold("(c) Scale Parameter"~b)),
     xlab = expression(b),
     ylab = "Frequency",
     col = col_hist2,
     border = "white",
     lwd = 0.5,
     xlim = range(c(bm_params[, 2], mean_scale_regular)) * c(0.95, 1.05)
)

# Add density overlay
d2 <- density(bm_params[, 2], bw = "SJ")
par(new = TRUE)
plot(d2$x, d2$y * diff(h2$mids[1:2]) * nrow(bm_params),
     type = "l",
     lwd = 2,
     col = "gray30",
     xlim = range(c(bm_params[, 2], mean_scale_regular)) * c(0.95, 1.05),
     ylim = c(0, max(h2$counts)),
     axes = FALSE,
     xlab = "",
     ylab = "")

# Add vertical lines for means
abline(v = mean_scale_bm, col = col_bm, lwd = 2, lty = 1)
abline(v = mean_scale_regular, col = col_regular, lwd = 2, lty = 2)

# Add mean values as text
text(mean_scale_bm, max(h2$counts) * 0.95,
     paste0("BM: ", round(mean_scale_bm, 4)),
     pos = 2, cex = 0.9, font = 2, col = col_bm)
text(mean_scale_regular, max(h2$counts) * 0.85,
     paste0("True: ", round(mean_scale_regular, 4)),
     pos = 4, cex = 0.9, font = 2, col = col_regular)

# box(lwd = 1.2, col = "gray40")

# ---- Panel D: Shape Parameter Distribution ----
if(ncol(bm_params) == 3) {
  h3 <- hist(bm_params[, 3], 
             breaks = 25,
             plot = FALSE)
  
  plot(h3,
       main = expression(bold("(c) Shape Parameter"~sigma)),
       xlab = expression(sigma),
       ylab = "Frequency",
       col = col_hist2,
       border = "white",
       lwd = 0.5,
       xlim = range(c(bm_params[, 3], mean_shape_regular)) * c(0.95, 1.05)
  )
  
  # Add density overlay
  d3 <- density(bm_params[, 2], bw = "SJ")
  par(new = TRUE)
  plot(d3$x, d3$y * diff(h3$mids[1:2]) * nrow(bm_params),
       type = "l",
       lwd = 2.5,
       col = "gray30",
       xlim = range(c(bm_params[, 3], mean_shape_regular)) * c(0.95, 1.05),
       ylim = c(0, max(h3$counts)),
       axes = FALSE,
       xlab = "",
       ylab = "")
  
  # Add vertical lines for means
  abline(v = mean_shape_bm, col = col_bm, lwd = 2, lty = 1)
  abline(v = mean_shape_regular, col = col_regular, lwd = 2, lty = 2)
  
  # Add mean values as text
  text(mean_shape_bm, max(h3$counts) * 0.95,
       paste0("BM: ", round(mean_shape_bm, 4)),
       pos = 2, cex = 0.9, font = 2, col = col_bm)
  text(mean_shape_regular, max(h3$counts) * 0.85,
       paste0("True: ", round(mean_shape_regular, 4)),
       pos = 4, cex = 0.9, font = 2, col = col_regular)
  
  box(lwd = 1.2, col = "gray40")
}

# ============================================
# Add summary statistics as a table below
# ============================================

# # Reset plot parameters for summary
# par(mfrow = c(1, 1), mar = c(0.5, 0.5, 2, 0.5))
# plot.new()
# 
# # Create summary text
# summary_text <- paste0(
#   "Summary Statistics (M = ", M, " replications, N = ", N, " observations, ", 
#   block_count, " blocks):\n",
#   "─────────────────────────────────────────────────────────────────────\n",
#   sprintf("%-25s %12s %12s %12s\n", "Method", "Location (μ)", "Scale (σ)", "Ratio"),
#   "─────────────────────────────────────────────────────────────────────\n",
#   sprintf("%-25s %12.4f %12.4f %12s\n", "Standard", mean_loc_regular, mean_scale_regular, "—"),
#   sprintf("%-25s %12.4f %12.4f %12.4f\n", "Block Maxima", mean_loc_bm, mean_scale_bm, mean(ratios[, 2])),
#   sprintf("%-25s %12.4f %12.4f %12.4f\n", "Block Maxima (adjusted)", 
#           mean(adjusted_locations), mean_scale_bm, mean(ratios[, 1])),
#   "─────────────────────────────────────────────────────────────────────\n"
# )
# 
# # Display summary with monospace font
# text(0.5, 0.7, summary_text, 
#      family = "mono", 
#      cex = 0.8, 
#      adj = c(0.5, 1))

# ============================================
# Save to file
# ============================================

# Save as PDF (recommended for publications)
dev.copy(pdf, "paper/plots/figA2", width = 8, height = 3)
dev.off()

# Save as high-resolution PNG
# dev.copy(png, "gumbel_comparison_publication.png",
#          width = 14, height = 5, units = "in", res = 300)
# dev.off()

# Save as EPS (often required by journals)
# dev.copy2eps(file = "gumbel_comparison_publication.eps", width = 14, height = 5)

cat("\n=== Publication Figure Generated ===\n")
cat("Figure dimensions: 14 x 5 inches\n")
cat("Panels: (a) Density comparison, (b) Location parameter, (c) Scale parameter\n")
cat("Color scheme: Colorblind-friendly palette\n")