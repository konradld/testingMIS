library(dplyr); library(evd)
source('R/4_bootstrap-dfb.R')

# ── Data & Model ──────────────────────────────────────────────────────────────
load("data/lsa.RData")

lsa <- lsa |> mutate(
  race1   = relevel(race1, ref="white"),
  fam_inc = relevel(factor(fam_inc), ref="3"),
  decile1 = relevel(factor(decile1), ref="5"),
  decile3 = relevel(factor(decile3), ref="5")
)

f <- as.integer(bar) ~ race1 + gender + factor(fam_inc) + factor(cluster) +
  lsat + ugpa + factor(decile1) + factor(decile3) + I(98 - age)

mdl <- lm(f, data=lsa)
summary(mdl)

# ── FWL Orthogonalisation ─────────────────────────────────────────────────────
mm      <- model.matrix(mdl)
y       <- mdl$model[, 1]
X1      <- mm[, 4]        # race1Black
Xother  <- mm[, -4]

fwl_vars <- fwl(y=y, X1=X1, X2=Xother)
fwl_lm   <- lm(fwl_vars[,1] ~ fwl_vars[,2] - 1)

block_count <- 30

# ── Analysis helper ───────────────────────────────────────────────────────────
analyse <- function(S) {
  res     <- estimate_dfb_evd(y, X1, Xother, S, block_count=block_count)
  p       <- res$params
  loc_adj <- p[1] + p[2] * log(block_count)
  Sdfb    <- abs(res$set_dfb)
  
  plot(density(res$block_maxima), col=0,
       xlim=c(0, max(Sdfb, max(res$block_maxima))), main="", xlab="")
  curve(dgumbel(x, loc_adj, p[2]), add=TRUE, col=4, lwd=2, lty=2)
  curve(dgumbel(x, p[1],    p[2]), add=TRUE, col=4, lwd=2, lty=2)
  abline(v=Sdfb)
  
  cat("p-value:", 1 - pgumbel(Sdfb, loc_adj, p[2]), "\n")
  
  fwl_lm2 <- lm(fwl_vars[-S,1] ~ fwl_vars[-S,2] - 1)
  cat("Coef shift:", fwl_lm2$coefficients - fwl_lm$coefficients, "\n")
  summary(fwl_lm2)
}

# ── Negative extremes ─────────────────────────────────────────────────────────
analyse(S = which(dfbeta(fwl_lm) < -.001))

summary(lm(f, data=lsa[-which(dfbeta(fwl_lm) < -.001), ]))

# ── Positive extremes ─────────────────────────────────────────────────────────
analyse(S = which(dfbeta(fwl_lm) > .0006))