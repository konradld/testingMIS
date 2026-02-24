library("testingMIS")
# devtools::load_all()
library("dplyr")

# ── Data & Model ──────────────────────────────────────────────────────────────
col_names <- c(
  "state",
  "county",
  "community",
  "communityname",
  "fold",
  "population",
  "householdsize",
  "racePctBlack",
  "racePctWhite",
  "racePctAsian",
  "racePctHisp",
  "agePct12t21",
  "agePct12t29",
  "agePct16t24",
  "agePct65up",
  "numbUrban",
  "pctUrban",
  "medIncome",
  "pctWWage",
  "pctWFarmSelf",
  "pctWInvInc",
  "pctWSocSec",
  "pctWPubAsst",
  "pctWRetire",
  "medFamInc",
  "perCapInc",
  "whitePerCap",
  "blackPerCap",
  "indianPerCap",
  "AsianPerCap",
  "OtherPerCap",
  "HispPerCap",
  "NumUnderPov",
  "PctPopUnderPov",
  "PctLess9thGrade",
  "PctNotHSGrad",
  "PctBSorMore",
  "PctUnemployed",
  "PctEmploy",
  "PctEmplManu",
  "PctEmplProfServ",
  "PctOccupManu",
  "PctOccupMgmtProf",
  "MalePctDivorce",
  "MalePctNevMarr",
  "FemalePctDiv",
  "TotalPctDiv",
  "PersPerFam",
  "PctFam2Par",
  "PctKids2Par",
  "PctYoungKids2Par",
  "PctTeen2Par",
  "PctWorkMomYoungKids",
  "PctWorkMom",
  "NumIlleg",
  "PctIlleg",
  "NumImmig",
  "PctImmigRecent",
  "PctImmigRec5",
  "PctImmigRec8",
  "PctImmigRec10",
  "PctRecentImmig",
  "PctRecImmig5",
  "PctRecImmig8",
  "PctRecImmig10",
  "PctSpeakEnglOnly",
  "PctNotSpeakEnglWell",
  "PctLargHouseFam",
  "PctLargHouseOccup",
  "PersPerOccupHous",
  "PersPerOwnOccHous",
  "PersPerRentOccHous",
  "PctPersOwnOccup",
  "PctPersDenseHous",
  "PctHousLess3BR",
  "MedNumBR",
  "HousVacant",
  "PctHousOccup",
  "PctHousOwnOcc",
  "PctVacantBoarded",
  "PctVacMore6Mos",
  "MedYrHousBuilt",
  "PctHousNoPhone",
  "PctWOFullPlumb",
  "OwnOccLowQuart",
  "OwnOccMedVal",
  "OwnOccHiQuart",
  "RentLowQ",
  "RentMedian",
  "RentHighQ",
  "MedRent",
  "MedRentPctHousInc",
  "MedOwnCostPctInc",
  "MedOwnCostPctIncNoMtg",
  "NumInShelters",
  "NumStreet",
  "PctForeignBorn",
  "PctBornSameState",
  "PctSameHouse85",
  "PctSameCity85",
  "PctSameState85",
  "LemasSwornFT",
  "LemasSwFTPerPop",
  "LemasSwFTFieldOps",
  "LemasSwFTFieldPerPop",
  "LemasTotalReq",
  "LemasTotReqPerPop",
  "PolicReqPerOffic",
  "PolicPerPop",
  "RacialMatchCommPol",
  "PctPolicWhite",
  "PctPolicBlack",
  "PctPolicHisp",
  "PctPolicAsian",
  "PctPolicMinor",
  "OfficAssgnDrugUnits",
  "NumKindsDrugsSeiz",
  "PolicAveOTWorked",
  "LandArea",
  "PopDens",
  "PctUsePubTrans",
  "PolicCars",
  "PolicOperBudg",
  "LemasPctPolicOnPatr",
  "LemasGangUnitDeploy",
  "LemasPctOfficDrugUn",
  "PolicBudgPerPop",
  "ViolentCrimesPerPop"
)

data <- read.csv("data/communities.data", header = FALSE, na.strings = "?") |>
  setNames(col_names)

clean_data <- \(df) {
  df[, colSums(is.na(df)) == 0] |>
    dplyr::select(-communityname, -state, -fold)
}

fit_mdl <- \(df) lm(ViolentCrimesPerPop ~ ., data = clean_data(df))

mdl <- fit_mdl(data)
mdl_rm <- fit_mdl(data[-c(463, 1583), ])
summary(mdl)
summary(mdl_rm)

sm <- summary(mdl)
(summary(mdl_rm)$coefficients[, 4] - sm$coefficients[, 4])[
  rownames(sm$coefficients[sm$coefficients[, 4] < 0.05, ])
]

# ── EVD Setup ─────────────────────────────────────────────────────────────────
mm <- model.matrix(mdl)
y <- as.matrix(mdl$model[, 1])
X1 <- mm[, "racePctBlack"]
Xother <- mm[, colnames(mm) != "racePctBlack"]

evd_plot <- function(y, X1, Xother, S, block_count, title) {
  res <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
  p <- res$params
  loc_adj <- p[1] + p[2] * log(block_count)
  Sdfb <- abs(res$set_dfb)

  plot(
    density(res$block_maxima),
    col = 0,
    xlim = c(0, max(Sdfb, max(res$block_maxima))),
    main = title,
    xlab = "",
    ylab = ""
  )
  curve(evd::dgumbel(x, loc_adj, p[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
  abline(v = Sdfb, col = 2, lwd = 2)
  cat(title, "p-value:", 1 - evd::pgumbel(Sdfb, loc_adj, p[2]), "\n")
}

# ── Analyses ──────────────────────────────────────────────────────────────────
par(mfrow = c(1, 3))
S1 <- c(463, 1583)
S2 <- c(1271, 669)

evd_plot(y, X1, Xother, c(S1, S2), block_count = 40, "Full Set")
evd_plot(y, X1, Xother, S1, block_count = 40, "First Set")

# Second set after removing first
evd_plot(
  y[-S1],
  X1[-S1],
  Xother[-S1, ],
  S = c(1271 - sum(S1 < 1271), 669 - sum(S1 < 669)), # adjust indices
  block_count = 40,
  "Second Set"
)

# ── Coefficient Comparison ────────────────────────────────────────────────────
mdl_rm2 <- fit_mdl(data[-c(S1, S2), ])
coefs <- sapply(list(mdl, mdl_rm, mdl_rm2), \(m) {
  summary(m)$coefficients["racePctBlack", 1]
}) |>
  setNames(c("orig", "part1", "all"))

cat(sprintf(
  "Shifts:  part1-orig = %.4f | all-orig = %.4f | all-part1 = %.4f\n",
  coefs["part1"] - coefs["orig"],
  coefs["all"] - coefs["orig"],
  coefs["all"] - coefs["part1"]
))
cat(sprintf(
  "Ratios:  all/orig = %.4f | part1/orig = %.4f | part1/all = %.4f\n",
  coefs["all"] / coefs["orig"],
  coefs["part1"] / coefs["orig"],
  coefs["part1"] / coefs["all"]
))
