library(evd)
source('R/01_estimate_dfb_evd.R')

# ── Data & Model ──────────────────────────────────────────────────────────────
sparrow <- read.delim("data/SparrowsElphick.txt")
f       <- tarsus ~ head + wingcrd + flatwing + culmen + nalospi + wt + Sex + Age
mdl     <- lm(f, data = sparrow)

headtail <- \(x) c(head=head(x), tail=tail(x))
dfbeta(mdl)[, 2] |> sort() |> headtail()
dfbeta(lm(f, data = sparrow[-646, ]))[, 2] |> sort() |> headtail()

# ── Setup ─────────────────────────────────────────────────────────────────────
y      <- mdl$model[, 1]
X1     <- model.matrix(mdl)[, "head"]
Xother <- model.matrix(mdl)[, colnames(model.matrix(mdl)) != "head"]

# ── Analysis helper ───────────────────────────────────────────────────────────
analyse <- function(S, block_count, shape = 0) {
  res    <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
  p      <- res$params
  loc_adj <- p[1] + p[2] * log(block_count)
  Sdfb   <- abs(res$set_dfb)
  
  plot(density(res$block_maxima), col=0,
       xlim=c(0, max(Sdfb, max(res$block_maxima))), main="", xlab="")
  curve(dgumbel(x, loc_adj, p[2]), add=TRUE, col=4, lwd=2, lty=2)
  abline(v=Sdfb)
  
  pval <- if (shape == 0) pgumbel(Sdfb, loc_adj, p[2])
  else            pgev(Sdfb,    loc_adj, p[2], shape=shape)
  cat("S =", S, "| p-value =", 1 - pval, "\n")
}

# ── Analyses ──────────────────────────────────────────────────────────────────
analyse(S=646, block_count=25)
analyse(S=366, block_count=40, shape=0.3183)
analyse(S=c(366, 646), block_count=25)