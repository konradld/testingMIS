library(evd)
source('R/01_estimate_dfb_evd.R')
source('R/06_fwl.R')

# ── Data & Model ──────────────────────────────────────────────────────────────
d   <- MASS::Boston
mdl <- lm(medv ~ ., data=d)
summary(mdl)

# ── FWL Orthogonalisation ─────────────────────────────────────────────────────
mm     <- model.matrix(mdl)
y      <- mdl$model[, 1]
X1     <- mm[, 2]        # crim
Xother <- mm[, -2]

fwl_vars <- fwl(y=y, X1=X1, X2=Xother)
fwl_lm   <- lm(fwl_vars[,1] ~ fwl_vars[,2] - 1)

dfb_ortho <- dfbeta(fwl_lm)
which(dfb_ortho %in% head(sort(dfb_ortho)))

# ── Analysis ──────────────────────────────────────────────────────────────────
S           <- c(381, 419, 406, 411, 366, 369)
block_count <- 35

res     <- estimate_dfb_evd(y, X1, Xother, S, block_count=block_count)
p       <- res$params
loc_adj <- p[1] + p[2] * log(block_count)
Sdfb    <- abs(res$set_dfb)

plot(density(res$block_maxima), col=0,
     xlim=c(0, max(Sdfb, max(res$block_maxima))), main="", xlab="")
curve(dgumbel(x, loc_adj, p[2]), add=TRUE, col=4, lwd=2, lty=2)
curve(dgumbel(x, p[1],    p[2]), add=TRUE, col=4, lwd=2, lty=2)
abline(v=Sdfb)

cat("p-value:", 1 - pgev(Sdfb, loc_adj, p[2], shape=0.2916), "\n")

fwl_lm2 <- lm(fwl_vars[-S,1] ~ fwl_vars[-S,2] - 1)
cat("Coef shift:", fwl_lm2$coefficients - fwl_lm$coefficients, "\n")