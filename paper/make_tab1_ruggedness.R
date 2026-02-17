# EVD Influence Analysis Summary Table Generation
# Load required packages
library(dplyr)
library(stargazer)
library(evd)

ruggedness <- read.csv("https://github.com/nk027/influential_sets/raw/refs/heads/main/paper/data/rugged_data.csv") |>
  dplyr::mutate(diamonds = gemstones / (land_area / 100))

# Prepare data
mf <- model.frame(log(rgdppc_2000) ~ rugged * cont_africa +
                    dist_coast * cont_africa, data = ruggedness)
y <- model.response(mf)
X <- model.matrix(mf, data = ruggedness)
mdl <- lm(y ~ X - 1)

source('R/4_bootstrap-dfb.R')

# Prepare orthogonalized variables
lm_data <- model.matrix(mdl)
y = mdl$model[,1]
X1 = lm_data[,5]  # ruggedness coefficient
Xother = lm_data[,-5]
fwl_vars <- fwl(y = y, X1 = X1, X2 = Xother)
yortho <- fwl_vars[,1]
Xortho <- fwl_vars[,2]
fwl_lm <- lm(yortho~Xortho-1)
block_count = 35

# Define country sets to test
country_sets <- list(
  "Seychelles" = c("Seychelles"),
  "Seychelles + Lesotho" = c("Seychelles", "Lesotho"),
  "Seychelles + Rwanda" = c("Seychelles", "Rwanda"),
  "Seychelles + Eswatini" = c("Seychelles", "Eswatini"),
  "Seychelles + Comoros" = c("Seychelles", "Comoros")#,
  # "Seychelles + Lesotho + Rwanda" = c("Seychelles", "Lesotho", "Rwanda"),
  # "Seychelles + Rwanda + Eswatini + Comoros" = c("Seychelles", "Lesotho","Eswatini", "Rwanda", "Comoros")
)

# Function to run EVD analysis for a country set
run_evd_analysis <- function(countries) {
  # Find row indices for these countries
  country_indices <- which(ruggedness[,3] %in% countries)
  S <- which(rownames(mf) %in% as.character(country_indices))
  
  # Run EVD analysis
  dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
  
  # Extract results
  params <- dfb_evd$params
  Sdfb <- dfb_evd$set_dfb
  
  # Calculate p-value
  p_value <- 1 - evd::pgumbel(abs(Sdfb), 
                              loc = params[1] + params[2] * log(block_count), 
                              scale = params[2])
  
  return(list(
    countries = paste(countries, collapse = " + "),
    mu = round(params[1], 4),
    sigma = round(params[2], 4),
    adj_location = round(params[1] + params[2] * log(block_count), 4),
    dfbeta = round(abs(Sdfb), 4),
    p_value = round(p_value, 4)
  ))
}

# Run analysis for each country set
results_list <- lapply(country_sets, run_evd_analysis)

# Create summary data frame
summary_df <- data.frame(
  # Set = paste0("Set ", 1:length(results_list)),
  Countries = sapply(results_list, function(x) x$countries),
  `|DFBETA|` = sapply(results_list, function(x) x$dfbeta),
  `Adjusted Location` = sapply(results_list, function(x) x$adj_location),
  `EVD Scale` = sapply(results_list, function(x) x$sigma),
  `P-value` = sapply(results_list, function(x) x$p_value),
  check.names = FALSE
)

# Print the table
print(summary_df)

# Create stargazer table
stargazer(summary_df, 
          type = "text",
          title = "EVD Influence Analysis Summary",
          summary = FALSE,
          rownames = FALSE,
          header = FALSE,
          notes = c("M = 35"))

# For LaTeX output
stargazer(summary_df, 
          type = "latex",
          title = "EVD Influence Analysis Summary",
          summary = FALSE,
          rownames = FALSE,
          header = FALSE,
          label = "tab:evd_influence",
          notes = c("M = 35"))
