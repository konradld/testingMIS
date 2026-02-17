rm(list = ls())
try(dev.off())

# set.seed(1234)
source('R/4_bootstrap-dfb.R')

# Load required libraries
library(evd)
library(extRemes)
library(stargazer)

# ============================================
# PARAMETERS
# ============================================
nS <- 1
N_vec <- c(20, 30, 50, 75, 100, 200, 300, 500, 1000) + nS
sx <- sqrt(2)
sr <- sqrt(2)
Rep <- 1000  # Number of replications for averaging

# ============================================
# FUNCTION DEFINITIONS
# ============================================

# Function to generate samples and fit GEV distribution
fit_theor_evd <- function(xdist, rdist, N, nS) {
  success <- FALSE
  while (!success) {
    fit_evd <- tryCatch({
      maxdraws <- replicate(N - nS, rmaxdfbeta(n = N - nS, xdist = xdist, 
                                               rdist = rdist, nS = nS))
      tmp_res <- extRemes::fevd(maxdraws)
      
      success <- TRUE  
      if(any(is.na(tmp_res$results$par))) { 
        success <- FALSE
        cat("NA in Coef...\n") 
      }
      tmp_res  
    }, error = function(e) {
      cat("Error occurred:", e$message, "\nRetrying...\n")
      NULL  
    })
  }
  param <- fit_evd$results$par
  return(param)
}

# Function to run simulation for all distributions at a given N
run_simulation_for_N <- function(N, nS, sx, sr, Rep) {
  cat("\n=== Running simulations for N =", N, "===\n")
  
  # Normal-Normal
  cat("  Normal-Normal...")
  xdist <- \(n) rnorm(n, 0, sx)
  rdist <- \(n) rnorm(n, 0, sr)
  thin_params <- replicate(Rep, fit_theor_evd(xdist, rdist, N, nS))
  cat(" Done\n")
  
  # t(5)-Normal
  cat("  t(5)-Normal...")
  xdist <- \(n) rt(n, 5)
  rdist <- \(n) rnorm(n, 0, sr)
  heavy_params1 <- replicate(Rep, fit_theor_evd(xdist, rdist, N, nS))
  cat(" Done\n")
  
  # Normal-t(5)
  cat("  Normal-t(5)...")
  xdist <- \(n) rnorm(n, 0, sr)
  rdist <- \(n) rt(n, 5)
  heavy_params2 <- replicate(Rep, fit_theor_evd(xdist, rdist, N, nS))
  cat(" Done\n")
  
  # t(5)-t(5)
  cat("  t(5)-t(5)...")
  xdist <- \(n) rt(n, 5)
  rdist <- \(n) rt(n, 5)
  fat_params <- replicate(Rep, fit_theor_evd(xdist, rdist, N, nS))
  cat(" Done\n")
  
  return(list(
    thin = thin_params,
    heavy1 = heavy_params1,
    heavy2 = heavy_params2,
    fat = fat_params
  ))
}

# Function to extract summary statistics
extract_summary_stats <- function(params_list, N) {
  dist_names <- c("Normal-Normal", "t(5)-Normal", "Normal-t(5)", "t(5)-t(5)")
  param_types <- c("thin", "heavy1", "heavy2", "fat")
  
  results <- data.frame()
  
  for (i in 1:4) {
    shape_vec <- params_list[[param_types[i]]][3, ]
    
    row_data <- data.frame(
      N = N,
      Distribution = dist_names[i],
      Mean_Shape = mean(shape_vec),
      SD_Shape = sd(shape_vec),
      Median_Shape = median(shape_vec),
      Q05_Shape = quantile(shape_vec, 0.05),
      Q25_Shape = quantile(shape_vec, 0.25),
      Q75_Shape = quantile(shape_vec, 0.75),
      Q95_Shape = quantile(shape_vec, 0.95)
    )
    
    results <- rbind(results, row_data)
  }
  
  return(results)
}

# ============================================
# MAIN SIMULATION LOOP
# ============================================

all_results <- data.frame()
all_params_storage <- list()  # Store all params for potential later use

for (N in N_vec) {
  params_list <- run_simulation_for_N(N, nS, sx, sr, Rep)
  summary_stats <- extract_summary_stats(params_list, N)
  all_results <- rbind(all_results, summary_stats)
  
  # Store full parameter results
  all_params_storage[[paste0("N_", N)]] <- params_list
}

# ============================================
# CREATE COMPREHENSIVE STARGAZER TABLE
# ============================================

# Reshape data for better table presentation
# Option 1: Wide format with N as columns

# Create separate tables for each statistic
table_mean <- reshape(all_results[, c("N", "Distribution", "Mean_Shape")],
                      idvar = "Distribution", timevar = "N", direction = "wide")
table_sd <- reshape(all_results[, c("N", "Distribution", "SD_Shape")],
                    idvar = "Distribution", timevar = "N", direction = "wide")
table_median <- reshape(all_results[, c("N", "Distribution", "Median_Shape")],
                        idvar = "Distribution", timevar = "N", direction = "wide")

# Clean column names
colnames(table_mean) <- c("Distribution", paste0("N=", N_vec))
colnames(table_sd) <- c("Distribution", paste0("N=", N_vec))
colnames(table_median) <- c("Distribution", paste0("N=", N_vec))

# Set row names
rownames(table_mean) <- table_mean$Distribution
rownames(table_sd) <- table_sd$Distribution
rownames(table_median) <- table_median$Distribution

# Remove Distribution column from data
table_mean <- table_mean[, -1]
table_sd <- table_sd[, -1]
table_median <- table_median[, -1]

# ============================================
# COMPREHENSIVE TABLE: All Statistics
# ============================================

# Create a comprehensive table with all key statistics for selected N values
# (e.g., N = 101, 301, 501, 701, 901, 1001)
selected_N <- c(101, 301, 501, 701, 901, 1001)
selected_results <- all_results[all_results$N %in% selected_N, ]

# Create formatted table
comp_table <- data.frame(
  Distribution = selected_results$Distribution,
  N = selected_results$N,
  Mean = selected_results$Mean_Shape,
  SD = selected_results$SD_Shape,
  Median = selected_results$Median_Shape,
  Q05 = selected_results$Q05_Shape,
  Q95 = selected_results$Q95_Shape
)

cat("\n=== COMPREHENSIVE SUMMARY TABLE (Selected N) ===\n")
stargazer(comp_table,
          type = "text",
          title = "GEV Shape Parameter Summary Statistics for Selected Sample Sizes",
          summary = FALSE,
          digits = 4,
          digit.separator = "",
          rownames = FALSE,
          notes = c(paste0("Based on ", Rep, " replications."),
                    "Q05 and Q95 are 5th and 95th percentiles."),
          notes.align = "l")

# LaTeX version
stargazer(comp_table,
          type = "latex",
          title = "GEV Shape Parameter Summary Statistics for Selected Sample Sizes",
          label = "tab:comprehensive_shape_params",
          summary = FALSE,
          digits = 4,
          digit.separator = "",
          rownames = FALSE,
          notes = c(paste0("Based on ", Rep, " replications."),
                    "$Q_{05}$ and $Q_{95}$ are 5th and 95th percentiles."),
          notes.align = "l",
          out = "paper/tables/comprehensive_shape_parameters_table.tex")

# ============================================
# FULL RESULTS TABLE
# ============================================

cat("\n=== FULL RESULTS TABLE (All N) ===\n")
full_table <- data.frame(
  N = all_results$N,
  Distribution = all_results$Distribution,
  Mean = all_results$Mean_Shape,
  SD = all_results$SD_Shape,
  Q25 = all_results$Q25_Shape,
  Median = all_results$Median_Shape,
  Q75 = all_results$Q75_Shape
)

stargazer(full_table,
          type = "text",
          title = "Complete GEV Shape Parameter Statistics across All Sample Sizes",
          summary = FALSE,
          digits = 4,
          digit.separator = "",
          rownames = FALSE,
          notes = paste0("Based on ", Rep, " replications for each N and distribution."),
          notes.align = "l")

# LaTeX version
stargazer(full_table,
          type = "latex",
          title = "Complete GEV Shape Parameter Statistics across All Sample Sizes",
          label = "tab:full_shape_params",
          summary = FALSE,
          digits = 4,
          digit.separator = "",
          rownames = FALSE,
          notes = paste0("Based on ", Rep, " replications for each $N$ and distribution."),
          notes.align = "l",
          out = "paper/tables/full_shape_parameters_table.tex")

# ============================================
# STATISTICAL TESTS
# ============================================

cat("\n=== Statistical Tests: H0: ξ⁻¹ = 0 (Gumbel) ===\n\n")

test_results <- data.frame()

for (N in N_vec) {
  params <- all_params_storage[[paste0("N_", N)]]
  
  # Extract shape parameters
  shape_thin <- params$thin[3, ]
  shape_heavy1 <- params$heavy1[3, ]
  shape_heavy2 <- params$heavy2[3, ]
  shape_fat <- params$fat[3, ]
  
  # T-tests
  test_thin <- t.test(shape_thin, mu = 0)
  test_heavy1 <- t.test(shape_heavy1, mu = 0)
  test_heavy2 <- t.test(shape_heavy2, mu = 0)
  test_fat <- t.test(shape_fat, mu = 0)
  
  # Store results
  test_results <- rbind(test_results, data.frame(
    N = N,
    Dist = "N-N",
    t_stat = test_thin$statistic,
    p_value = test_thin$p.value
  ))
  test_results <- rbind(test_results, data.frame(
    N = N,
    Dist = "t-N",
    t_stat = test_heavy1$statistic,
    p_value = test_heavy1$p.value
  ))
  test_results <- rbind(test_results, data.frame(
    N = N,
    Dist = "N-t",
    t_stat = test_heavy2$statistic,
    p_value = test_heavy2$p.value
  ))
  test_results <- rbind(test_results, data.frame(
    N = N,
    Dist = "t-t",
    t_stat = test_fat$statistic,
    p_value = test_fat$p.value
  ))
}

print(test_results)

# Save test results
# stargazer(test_results,
#           type = "latex",
#           title = "T-test Results for Shape Parameters (H0: $\\xi^{-1} = 0$)",
#           label = "tab:ttests",
#           summary = FALSE,
#           digits = 3,
#           digit.separator = "",
#           rownames = FALSE,
#           out = "shape_parameter_ttests.tex")

# ============================================
# SAVE RESULTS
# ============================================

# Save all results to RData file for later use
save(all_results, all_params_storage, test_results, N_vec, Rep,
     file = "paper/tables/tabA1.RData")

# ============================================
# CREATE PLOTS
# ============================================

create_convergence_plot <- function(all_results, 
                                    output_file = "paper/plots/fig3.pdf",
                                    width = 9, 
                                    height = 5.5) {
  
  # Extract data in the required format
  data <- data.frame(
    N = all_results$N,
    Distribution = all_results$Distribution,
    Mean = all_results$Mean_Shape,
    StdDev = all_results$SD_Shape
  )
  
  # Rename distributions to match your style
  data$Distribution <- gsub("Normal-Normal", "Normal--Normal", data$Distribution)
  data$Distribution <- gsub("t\\(5\\)-Normal", "t(5)--Normal", data$Distribution)
  data$Distribution <- gsub("Normal-t\\(5\\)", "Normal--t(5)", data$Distribution)
  data$Distribution <- gsub("t\\(5\\)-t\\(5\\)", "t(5)--t(5)", data$Distribution)
  
  # Calculate confidence interval bounds
  data$Lower <- data$Mean - data$StdDev
  data$Upper <- data$Mean + data$StdDev
  
  # Get unique N values and distributions
  N_vals <- unique(data$N)
  distributions <- c("Normal--Normal", "t(5)--Normal", "Normal--t(5)", "t(5)--t(5)")
  
  # Colorblind-friendly colors (matching your palette)
  colors <- c("#440154", "#31688E", "#35B779", "#D55E00")
  
  # Set up the plot with LOG X-AXIS
  plot(NULL, xlim = c(min(N_vals) * 0.85, max(N_vals) * 1.15), 
       ylim = c(-0.05, 0.35), 
       xlab = "Sample Size (N)", 
       ylab = "", 
       axes = FALSE, 
       log = "x")
  
  mtext(expression(~xi^-1), side = 2, line = 3, las = 1)
  
  # Add horizontal reference lines
  lines(x = c(min(N_vals) * 0.75, max(N_vals) * 1.25), 
        y = c(0, 0), lty = 2, lwd = 2, col = "black")
  lines(x = c(min(N_vals) * 0.75, max(N_vals) * 1.25), 
        y = c(0.2, 0.2), lty = 2, lwd = 2, col = "black")
  
  # Plot confidence intervals and lines for each distribution
  for (i in 1:4) {
    dist_data <- data[data$Distribution == distributions[i], ]
    
    # Create polygon for confidence interval
    x_poly <- c(dist_data$N, rev(dist_data$N))
    y_poly <- c(dist_data$Lower, rev(dist_data$Upper))
    
    # Add confidence interval as polygon
    polygon(x_poly, y_poly, 
            col = adjustcolor(colors[i], alpha = 0.2), 
            border = NA)
    
    # Add line
    lines(dist_data$N, dist_data$Mean, col = colors[i], lwd = 2)
    
    # Add points
    points(dist_data$N, dist_data$Mean, col = colors[i], pch = 16, cex = 1.2)
  }
  
  # Add axes
  axis(1)
  axis(2, las = 1)
  
  # Add legend
  legend("top", 
         legend = c("N–N", "t(5)–N", "N–t(5)", "t(5)–t(5)"),
         col = colors, 
         lwd = 2, 
         pch = 16, 
         cex = 0.90, 
         bg = NA, 
         box.col = NA,
         horiz = TRUE)
  
  # Save plot
  dev.copy(pdf, output_file, width = width, height = height)
  dev.off()
  
  cat("\nConvergence plot saved to:", output_file, "\n")
  
  # Return the formatted data for inspection
  invisible(data)
}

# ============================================
# Usage example (if run standalone)
# ============================================

if (exists("all_results")) {
  # Create the plot from simulation results
  plot_data <- create_convergence_plot(all_results)
} else {
  cat("Note: Load gev_analysis_results.RData first to create plots automatically\n")
  cat("Or run generalized_gev_analysis.R to generate results\n")
}

# ============================================
# Function to print summary table in text format
# ============================================

print_convergence_table <- function(all_results) {
  
  # Extract and format data
  summary_table <- all_results[, c("N", "Distribution", "Mean_Shape", "SD_Shape")]
  
  # Rename for clarity
  names(summary_table) <- c("N", "Distribution", "Mean(ξ⁻¹)", "SD(ξ⁻¹)")
  
  # Round to 4 decimal places
  summary_table[, 3:4] <- round(summary_table[, 3:4], 4)
  
  cat("\n" , paste(rep("=", 70), collapse = ""), "\n")
  cat("SHAPE PARAMETER CONVERGENCE SUMMARY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Print by distribution
  distributions <- unique(summary_table$Distribution)
  
  for (dist in distributions) {
    cat("\n", dist, ":\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    dist_data <- summary_table[summary_table$Distribution == dist, ]
    
    # Print table header
    cat(sprintf("%-10s %15s %15s\n", "N", "Mean(ξ⁻¹)", "SD(ξ⁻¹)"))
    cat(paste(rep("-", 45), collapse = ""), "\n")
    
    # Print rows
    for (i in 1:nrow(dist_data)) {
      cat(sprintf("%-10d %15.4f %15.4f\n", 
                  dist_data$N[i], 
                  dist_data$`Mean(ξ⁻¹)`[i], 
                  dist_data$`SD(ξ⁻¹)`[i]))
    }
  }
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")
  
  invisible(summary_table)
}

# Print table if results exist
if (exists("all_results")) {
  print_convergence_table(all_results)
}
