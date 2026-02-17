rm(list = ls())
headtail <- \(x) c(head = head(x), tail = tail(x))

d <- MASS::Boston

mdl <- lm(medv ~ ., data = d)
rm <- c(381, 419, 406, 411, 366)[1L] # First is crass, 2, 3, and 4 have decent leverage, 5 is crass
mdl2 <- lm(medv ~ ., data = d[-rm, ])

summary(mdl)
summary(mdl2)
#===== Estimation of EVD =========

# ---- looking at crime
source('R/4_bootstrap-dfb.R')
lm_data <- model.matrix(mdl)
lm_data |> dimnames() |> _[[2]][2]

y = mdl$model[,1]
X1 = lm_data[,2]
Xother = lm_data[,-2]

fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
yortho <- fwl_vars[,1]
Xortho <- fwl_vars[,2]

fwl_lm <- lm(yortho~Xortho-1)
dfb_ortho <- stats::dfbeta(fwl_lm)
which(dfb_ortho %in% head(sort(dfb_ortho)))

block_count = 35
par(mfrow = c(1,1))

S <- c(381, 419, 406, 411, 366, 369)

# S <- c(419, 406, 411)

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
curve(evd::dgumbel(x,loc = params[1] , scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgev(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2], shape = 0.2916 )

#dfb
fwl_lm2 <- lm(yortho[-S]~Xortho[-S]-1)
fwl_lm2$coefficients - fwl_lm$coefficients

fwl_lm |> summary()
fwl_lm2 |> summary()


mdl2 <- lm(medv ~ ., data = d[-S, ])

summary(mdl)
summary(mdl2)
