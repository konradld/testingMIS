library("testingMIS")
# devtools::load_all()
library("dplyr")

# ── Data & Model ──────────────────────────────────────────────────────────────
ruggedness <- read.csv("data/rugged_data.csv") |>
  dplyr::mutate(diamonds = gemstones / (land_area / 100))
# Replace old name in the dataset
ruggedness[ruggedness[, "country"] == "Swaziland", "country"] <- "Eswatini"

mf <- model.frame(
  log(rgdppc_2000) ~ rugged * cont_africa + dist_coast * cont_africa,
  data = ruggedness
)

# Create model frame AND matrix together
mdl <- lm(mf)

y <- model.response(model.frame(mdl))
X <- model.matrix(mdl)

X1 <- X[, 5]
Xother <- X[, -5]

# Orthogonalize
fwl_vars <- testingMIS:::fwl(y = y, X = X1, Z = Xother)
fwl_lm <- lm(fwl_vars[, 1] ~ fwl_vars[, 2] - 1)

# ── Analysis ──────────────────────────────────────────────────────────────────
block_count <- 35

country_sets <- list(
  "Seychelles" = c("Seychelles"),
  "Seychelles + Lesotho" = c("Seychelles", "Lesotho"),
  "Seychelles + Rwanda" = c("Seychelles", "Rwanda"),
  "Seychelles + Eswatini" = c("Seychelles", "Eswatini"),
  "Seychelles + Comoros" = c("Seychelles", "Comoros")
)

run_evd <- function(countries) {
  S <- which(
    rownames(mf) %in% as.character(which(ruggedness[, 3] %in% countries))
  )
  res <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
  p <- res$params
  loc_adj <- p[1] + p[2] * log(block_count)
  Sdfb <- abs(res$set_dfb)
  data.frame(
    `|DFBETA|` = round(Sdfb, 4),
    `Adjusted Location` = round(loc_adj, 4),
    `EVD Scale` = round(p[2], 4),
    `P-value` = round(1 - evd::pgumbel(Sdfb, loc_adj, p[2]), 4),
    check.names = FALSE
  )
}

summary_df <- cbind(
  Countries = names(country_sets),
  do.call(rbind, lapply(country_sets, run_evd))
)
print(summary_df)
