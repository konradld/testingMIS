rm(list = ls())
library(evd)
source('R/4_bootstrap-dfb.R')

# ── Data & Model ──────────────────────────────────────────────────────────────
data <- read.csv("data/adult.data", header=FALSE, na.strings="?") |>
  setNames(c("age","workclass","fnlwgt","education","education_num",
             "marital_status","occupation","relationship","race","sex",
             "capital_gain","capital_loss","hours_per_week","native_country","income"))
data$y <- as.integer(data$income == " >50K")

mdl <- lm(y ~ age + education_num + hours_per_week + capital_gain + capital_loss +
            workclass + marital_status + occupation + relationship +
            race + sex + native_country, data=data)

# ── FWL Orthogonalisation ─────────────────────────────────────────────────────
mm     <- model.matrix(mdl)
mm     <- mm[, colnames(mm) != "occupation Transport-moving"]  # drop perfect multicollinearity
y      <- mdl$model[, 1]
X1     <- mm[, 43]
Xother <- mm[, -43]

fwl_vars <- fwl(y=y, X1=X1, X2=Xother)
fwl_lm   <- lm(fwl_vars[,1] ~ fwl_vars[,2] - 1)

block_count <- 35

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
}

dfb <- sort(dfbeta(fwl_lm)[, 1])

par(mfrow=c(1, 2))
analyse(as.numeric(names(head(dfb, 325))))  # negative extremes
analyse(as.numeric(names(tail(dfb, 325))))  # positive extremes