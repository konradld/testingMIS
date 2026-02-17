rm(list = ls())
# Law School and Communities and Crime are in another script
headtail <- \(x) c(head = head(x), tail = tail(x))
# Adult ---
data <- read.csv(
  "data/adult.data",
  header = FALSE,
  na.strings = "?"
)
colnames(data) <- c(
  "age",
  "workclass",
  "fnlwgt",
  "education",
  "education_num",
  "marital_status",
  "occupation",
  "relationship",
  "race",
  "sex",
  "capital_gain",
  "capital_loss",
  "hours_per_week",
  "native_country",
  "income"
)
data$y <- ifelse(data$income == " >50K", 1, 0)

mdl <- lm(
  y ~
    age +
    education_num +
    hours_per_week +
    capital_gain +
    capital_loss +
    workclass +
    marital_status +
    occupation +
    relationship +
    race +
    sex +
    native_country,
  data = data
)


#===== Estimation of EVD =========
source('R/4_bootstrap-dfb.R')
lm_data <- model.matrix(mdl)[,-which(colnames(model.matrix(mdl)) == "occupation Transport-moving")] # exclude 34 bc of perfect multicolinearity
lm_data |> dimnames() |> _[[2]][43]

y = mdl$model[,1]
X1 = lm_data[,43]
Xother = lm_data[,-43]

fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
yortho <- fwl_vars[,1]
Xortho <- fwl_vars[,2]

fwl_lm <- lm(yortho~Xortho-1)

block_count = 35
par(mfrow = c(1,2))
#---- negative extremes
# S <- stats::dfbeta(mdl)[, 44] |> sort() |> head(325) |> names() |> as.numeric()
S <- stats::dfbeta(fwl_lm)[,1] |> sort() |> head(325) |> names() |> as.numeric()

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

#---- positive extremes
# S <- stats::dfbeta(mdl)[, 44] |> sort() |> tail(325) |> names() |> as.numeric()
S <- stats::dfbeta(fwl_lm)[,1] |> sort() |> tail(325) |> names() |> as.numeric()


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
