library(evd)
source('R/01_estimate_dfb_evd.R')
source('R/03_dfb_bmx.R')
source('R/05_rmaxdfbeta.R')

# ── Parameters ────────────────────────────────────────────────────────────────
nS          <- 1
N           <- 5000 + nS
block_count <- 30
M           <- 100
xdist       <- \(n) rnorm(n, 0, sqrt(2))
rdist       <- \(n) rnorm(n, 0, sqrt(2))

pal <- list(reg = "#0072B2", bm = "#009E73", bma = "#D55E00",
            h1 = "#A6D8F4", h2 = "#F0E442")

# ── Simulation ────────────────────────────────────────────────────────────────
fit_both <- function() {
  repeat {
    res <- tryCatch({
      # Method 1: theoretical maxima
      p1 <- fgev(replicate(1000, rmaxdfbeta(N-nS, xdist, rdist, nS)),
                 std.err=FALSE, warn.inf=FALSE, shape=0)$estimate
      
      # Method 2: block maxima
      X  <- xdist(N); Y <- 2*X + rdist(N)
      R  <- resid(lm(Y ~ X - 1))
      p2 <- fgev(abs(dfb_bmx(X, R, 1:nS, block_count)),
                 std.err=FALSE, warn.inf=FALSE, shape=0)$estimate
      
      list(regular=p1, bm=p2)
    }, error=\(e) { cat("Error:", e$message, "\nRetrying...\n"); NULL })
    if (!is.null(res)) return(res)
  }
}

cat("Fitting Gumbel distributions for", M, "replications...\n")
results <- replicate(M, fit_both(), simplify=FALSE)

reg <- do.call(rbind, lapply(results, `[[`, "regular"))
bm  <- do.call(rbind, lapply(results, `[[`, "bm"))
bma <- bm[,1] + bm[,2] * log(block_count)  # adjusted locations

# ── Average densities ─────────────────────────────────────────────────────────
x_seq <- seq(0, 0.01, length.out=1000)

avg_dens <- function(locs, scales) {
  rowMeans(mapply(\(l,s) dgumbel(x_seq, l, s), locs, scales))
}

d_reg <- avg_dens(reg[,1], reg[,2])
d_bm  <- avg_dens(bm[,1],  bm[,2])
d_bma <- avg_dens(bma,      bm[,2])

# ── Plotting helpers ──────────────────────────────────────────────────────────
par(mfrow=c(1, 3 + (ncol(bm)==3)), mar=c(4.5,4.5,3,2), mgp=c(2.8,0.8,0),
    las=1, cex.lab=1.15, cex.axis=0.95, cex.main=1.2)

hist_panel <- function(vals, ref_val, col_fill, col_v1, col_v2,
                       main, xlab, lab1="BM", lab2="True") {
  h <- hist(vals, breaks=25, plot=FALSE)
  xlim <- range(c(vals, ref_val)) * c(0.95, 1.05)
  plot(h, main=main, xlab=xlab, ylab="Frequency",
       col=col_fill, border="white", xlim=xlim)
  d <- density(vals, bw="SJ")
  par(new=TRUE)
  plot(d$x, d$y * diff(h$mids[1:2]) * length(vals), type="l", lwd=2,
       col="gray30", xlim=xlim, ylim=c(0, max(h$counts)), axes=FALSE, xlab="", ylab="")
  abline(v=c(mean(vals), ref_val), col=c(col_v1, col_v2), lwd=2, lty=1:2)
  text(mean(vals),  max(h$counts)*0.95, sprintf("%s: %.4f", lab1, mean(vals)),  pos=2, cex=0.9, font=2, col=col_v1)
  text(ref_val,     max(h$counts)*0.85, sprintf("%s: %.4f", lab2, ref_val),     pos=4, cex=0.9, font=2, col=col_v2)
}

# ── Panel A: Density comparison ───────────────────────────────────────────────
plot(x_seq, d_reg, type="l", lwd=2, col=pal$reg,
     main=expression(bold("(a) Average Gumbel Densities")),
     xlab=expression(Delta^max), ylab="Density",
     xlim=c(0, 0.003), ylim=c(0, max(d_reg, d_bm, d_bma)*1.05), axes=FALSE)
axis(1, at=c(0, 0.001, 0.002, 0.003)); axis(2, las=1)
lines(x_seq, d_bm,  col=pal$bm,  lwd=2, lty=2)
lines(x_seq, d_bma, col=pal$bma, lwd=2, lty=4)
legend("topright", legend=c("  True","  BM","  BM (adj.)"),
       col=unlist(pal[c("reg","bm","bma")]), lwd=2, lty=c(1,2,4),
       bty="n", cex=1, y.intersp=0.9, inset=0.1, seg.len=4)

# ── Panel B: Location ─────────────────────────────────────────────────────────
hist_panel(bma, mean(reg[,1]), pal$h1, pal$bma, pal$reg,
           expression(bold("(b) Location Parameter"~tilde(a))),
           expression(tilde(a)))

# ── Panel C: Scale ────────────────────────────────────────────────────────────
hist_panel(bm[,2], mean(reg[,2]), pal$h2, pal$bm, pal$reg,
           expression(bold("(c) Scale Parameter"~b)),
           expression(b))

# ── Panel D: Shape (if GEV) ───────────────────────────────────────────────────
if (ncol(bm) == 3)
  hist_panel(bm[,3], mean(reg[,3]), pal$h2, pal$bm, pal$reg,
             expression(bold("(d) Shape Parameter"~sigma)),
             expression(sigma))

dev.copy(pdf, "paper/plots/figA1.pdf", width=8, height=3); dev.off()