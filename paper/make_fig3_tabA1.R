rm(list = ls())
try(dev.off())

library(evd); library(extRemes); library(stargazer)
source('R/4_bootstrap-dfb.R')

# ── Parameters ────────────────────────────────────────────────────────────────
nS    <- 1
N_vec <- c(20, 30, 50, 75, 100, 200, 300, 500, 1000) + nS
Rep   <- 1000

dists <- list(
  "Normal--Normal" = list(x = \(n) rnorm(n, 0, sqrt(2)), r = \(n) rnorm(n, 0, sqrt(2))),
  "t(5)--Normal"   = list(x = \(n) rt(n, 5),             r = \(n) rnorm(n, 0, sqrt(2))),
  "Normal--t(5)"   = list(x = \(n) rnorm(n, 0, sqrt(2)), r = \(n) rt(n, 5)),
  "t(5)--t(5)"     = list(x = \(n) rt(n, 5),             r = \(n) rt(n, 5))
)

# ── EVD Fitting ───────────────────────────────────────────────────────────────
fit_theor_evd <- function(xdist, rdist, N, nS) {
  repeat {
    fit <- tryCatch({
      draws <- replicate(N - nS, rmaxdfbeta(N - nS, xdist, rdist, nS))
      res   <- extRemes::fevd(draws)
      if (anyNA(res$results$par)) { cat("NA in Coef...\n"); NULL } else res
    }, error = \(e) { cat("Error:", e$message, "\nRetrying...\n"); NULL })
    if (!is.null(fit)) return(fit$results$par)
  }
}

# ── Main Simulation ───────────────────────────────────────────────────────────
all_params <- lapply(setNames(N_vec, paste0("N_", N_vec)), \(N) {
  cat("\n=== N =", N, "===\n")
  lapply(dists, \(d) { cat(" ", deparse(substitute(d)), "..."); 
    res <- replicate(Rep, fit_theor_evd(d$x, d$r, N, nS))
    cat(" Done\n"); res })
})

# ── Summary Statistics ────────────────────────────────────────────────────────
all_results <- do.call(rbind, lapply(N_vec, \(N) {
  do.call(rbind, lapply(names(dists), \(dn) {
    s <- all_params[[paste0("N_", N)]][[dn]][3, ]
    data.frame(N, Distribution = dn,
               Mean = mean(s), SD = sd(s), Median = median(s),
               Q05 = quantile(s, .05), Q25 = quantile(s, .25),
               Q75 = quantile(s, .75), Q95 = quantile(s, .95),
               row.names = NULL)
  }))
}))

# ── Tables ────────────────────────────────────────────────────────────────────
make_table <- function(df, type, title, label = NULL, file = NULL) {
  stargazer(df, type = type, title = title, label = label,
            summary = FALSE, digits = 4, digit.separator = "",
            rownames = FALSE,
            notes = c(sprintf("Based on %d replications.", Rep),
                      if (type == "latex") "$Q_{05}$/$Q_{95}$: 5th/95th percentiles."
                      else "Q05/Q95: 5th/95th percentiles."),
            notes.align = "l", out = file)
}

tab_full <- all_results[, c("N","Distribution","Mean","SD","Q25","Median","Q75")]

make_table(tab_full, "text",  "Complete GEV Shape Parameter Statistics")
make_table(tab_full, "latex", "Complete GEV Shape Parameter Statistics",
           "tab:full_shape_params", "paper/tables/full_shape_parameters_table.tex")

# ── Save ──────────────────────────────────────────────────────────────────────
save(all_results, all_params, N_vec, Rep,
     file = "paper/tables/tabA1.RData")

# ── Convergence Plot ──────────────────────────────────────────────────────────
plot_data <- transform(all_results,
                       Lower = Mean - SD, Upper = Mean + SD,
                       Distribution = gsub("-", "-", gsub("\\(", "(", Distribution)))

colors <- c("#440154","#31688E","#35B779","#D55E00")
N_vals <- unique(plot_data$N)

plot(NULL, xlim = range(N_vals) * c(0.85, 1.15), ylim = c(-0.05, 0.35),
     xlab = "Sample Size (N)", ylab = "", axes = FALSE, log = "x")
mtext(expression(~xi^-1), side = 2, line = 3, las = 1)
abline(h = c(0, 0.2), lty = 2, lwd = 2)

for (i in seq_along(dists)) {
  d <- plot_data[plot_data$Distribution == names(dists)[i], ]
  polygon(c(d$N, rev(d$N)), c(d$Lower, rev(d$Upper)),
          col = adjustcolor(colors[i], .2), border = NA)
  lines(d$N, d$Mean, col = colors[i], lwd = 2)
  points(d$N, d$Mean, col = colors[i], pch = 16, cex = 1.2)
}

axis(1); axis(2, las = 1)
legend("top", legend = gsub("Normal", "N", names(dists)),
       col = colors, lwd = 2, pch = 16, cex = .9,
       bg = NA, box.col = NA, horiz = TRUE)

dev.copy(pdf, "paper/plots/fig3.pdf", width = 9, height = 5.5); dev.off()
