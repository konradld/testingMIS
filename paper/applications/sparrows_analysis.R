rm(list = ls())
headtail <- \(x) c(head = head(x), tail = tail(x))

# Sparrows ---
# Multiple reg of tarsus length ~ head (susceptible), wing, beak length
sparrow <- read.delim("data/SparrowsElphick.txt")

# f <- tarsus ~ head + wingcrd + culmen
# mdl <- lm(f, data = sparrow)
# summary(mdl)
# stats::dfbeta(mdl)[, 2] |> sort() |> headtail()
# mdl_i <- lm(f, data = sparrow[-646, ])
# summary(mdl_i)
# stats::dfbeta(mdl_i)[, 2] |> sort() |> headtail()

f <- tarsus ~ head + wingcrd + flatwing + culmen + nalospi + wt + Sex + Age
mdl <- lm(f, data = sparrow)
summary(mdl)
stats::dfbeta(mdl)[, 2] |> sort() |> headtail()
mdl_i <- lm(f, data = sparrow[-646, ])
summary(mdl_i)
stats::dfbeta(mdl_i)[, 2] |> sort() |> headtail()

#########################################################
# ---- Biggest Outlier
source('R/4_bootstrap-dfb.R')
lm_data <- sparrow[,c('tarsus', 'head', 'wingcrd', 'flatwing', 'culmen', 'nalospi', 'wt', 'Sex', 'Age')]
lm_data <- cbind(lm_data, 1)

lm_data <- model.matrix(mdl)

y = mdl$model[,1]
X1 = lm_data[,'head']
Xother = lm_data[,-which(colnames(lm_data) == "head")]

S  <- c(646) # select the set

block_count = 25

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])

#########################################################
# ---- Set of 
X1 |> sort() |> headtail()
S  <- c(366) # select the set

block_count = 40
dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
# dfb_evd <- estimate_dfb_evd(y[-646], X1[-646], Xother[-646,], S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgev(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2], shape = 0.3183)
# 1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])

#########################################################
# ---- Set of 2
X1 |> sort() |> headtail()
S  <- c(366, 646) # select the set

block_count = 25

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "", xlab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 1)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])
