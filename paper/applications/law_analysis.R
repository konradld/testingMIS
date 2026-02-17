rm(list = ls())
# Law School
# data("lsa", package = "qeML")
load("05_input/lsa.RData") # I was being lazy
headtail <- \(x) c(head = head(x), tail = tail(x))

lsa <- lsa |>
  dplyr::mutate(
    race1 = relevel(race1, ref = "white"),
    fam_inc = relevel(factor(fam_inc), ref = "3"),
    decile1 = relevel(factor(decile1), ref = "5"),
    decile3 = relevel(factor(decile3), ref = "5")
  )

# Q: Passing the bar vs features, where sex and race are protected
# Student's grade decile in year 1 and 3
# LSAT score
# Undergraduate, 1st law-school, cumulative law-school GPA
mdl <- lm(
  as.integer(bar) ~
    race1 +
    gender +
    factor(fam_inc) +
    factor(cluster) + # Law school quality
    lsat +
    ugpa +
    factor(decile1) + # 1st year GPA decile
    factor(decile3) +
    I(98 - age), # Age is coded as the birthyear, data is from 1998
  data = lsa
)
summary(mdl)
dfb <- stats::dfbeta(mdl)[, 2:6]

#===== Estimation of EVD =========
source('00_code/dfbeta_funcs.R')
lm_data <- model.matrix(mdl)
lm_data |> dimnames() |> _[[2]][4]

y = mdl$model[,1]
X1 = lm_data[,4]
Xother = lm_data[,-4]

fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
yortho <- fwl_vars[,1]
Xortho <- fwl_vars[,2]

fwl_lm <- lm(yortho~Xortho-1)

block_count = 30
par(mfrow = c(1,1))
#---- negative extremes
# S <- which(dfb[, 4] < -.0005)
S <- which(stats::dfbeta(fwl_lm) < -.001)

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
curve(evd::dgumbel(x,loc = params[1] , scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])

#dfb
fwl_lm2 <- lm(yortho[-S]~Xortho[-S]-1)
fwl_lm2$coefficients - fwl_lm$coefficients

mdl2 <- lm(
  as.integer(bar) ~
    race1 +
    gender +
    factor(fam_inc) +
    factor(cluster) + # Law school quality
    lsat +
    ugpa +
    factor(decile1) + # 1st year GPA decile
    factor(decile3) +
    I(98 - age), # Age is coded as the birthyear, data is from 1998
  data = lsa[-S,]
)


fwl_lm |> summary()
fwl_lm2 |> summary()

mdl |> summary()
mdl2 |> summary()

#---- positive extremes
# S <- which(dfb[, 4] > .0005)
S <- which(stats::dfbeta(fwl_lm) > .0006)


dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
curve(evd::dgumbel(x,loc = params[1] , scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])

#dfb
fwl_lm2 <- lm(yortho[-S]~Xortho[-S]-1)
fwl_lm2$coefficients - fwl_lm$coefficients

fwl_lm |> summary()
fwl_lm2 |> summary()
