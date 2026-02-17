rm(list = ls())
data <- read.csv(
  "data/communities.data",
  header = FALSE,
  na.strings = "?"
)
colnames(data) <- c(
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

mdl <- lm(
  ViolentCrimesPerPop ~ .,
  data[, which(colSums(is.na(data)) == 0)] |>
    dplyr::select(-communityname, -state, -fold)
)
summary(mdl)
sm <- summary(mdl)

# Two & two arguably lower the "racePctBlack" coefficent
mdl_rm <- lm(
  ViolentCrimesPerPop ~ .,
  data[-c(463, 1583), which(colSums(is.na(data)) == 0)] |>
    dplyr::select(-communityname, -state, -fold)
)
summary(mdl_rm)

(summary(mdl_rm)$coefficients[,4] - summary(mdl)$coefficients[,4])[rownames(sm$coefficients[sm$coefficients[,4]<0.05,])]
# # Two and ... raise it by a bit
# mdl_rm <- lm(
#   ViolentCrimesPerPop ~ .,
#   data[-c(1271, 669), which(colSums(is.na(data)) == 0)] |>
#     dplyr::select(-communityname, -state, -fold)
# )
# summary(mdl_rm)

#########################################################
source('R/4_bootstrap-dfb.R')
# lm_data <- data[, which(colSums(is.na(data)) == 0)] |>
#   dplyr::select(-communityname, -state, -fold)
# lm_data <- cbind(lm_data, 1)

lm_data <- model.matrix(mdl)

y = mdl$model[1] |> as.matrix()
X1 = lm_data[,'racePctBlack']
Xother = lm_data[,!colnames(lm_data) %in% c('racePctBlack')]

# ---- Full Set
# S  <- c(463, 1583) # select the set
# S  <- c(1271, 669)
S <- c(463, 1583, 1271, 669)
block_count = 40

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

par(mfrow = c(1,2))
plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "First Set", xlab = "", ylab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 2, lwd = 2)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])

# ---- 1st Partial Set
S <- c(463, 1583)
block_count = 40

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

par(mfrow = c(1,2))
plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "First Set", xlab = "", ylab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 2, lwd = 2)

#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])



# ---- 2nd Partial Set
# after exclusion of the first set
y = y[-S]
X1 = X1[-S]
Xother = Xother[-S,]

# S  <- c(463-sum(S<463), 1583-sum(S<1583)) # select the set
S  <- c(1271-sum(S<1271), 669-sum(S<669))

block_count = 40

dfb_evd <- estimate_dfb_evd(y, X1, Xother, S, block_count = block_count)
params <- dfb_evd$params
Sdfb <- dfb_evd$set_dfb
D_bsmx <- dfb_evd$block_maxima

plot(density(D_bsmx), col = 0, xlim = c(0,max(abs(Sdfb),max(D_bsmx))), main = "Second Set", xlab = "", ylab = "")
curve(evd::dgumbel(x,loc = params[1] + params[2] * log(block_count), scale = params[2]), add = TRUE, col = 4, lwd = 2, lty = 2)
abline(v = abs(Sdfb), col = 2, lwd = 2)


#p-value
1-evd::pgumbel(abs(Sdfb),loc = params[1] + params[2] * log(block_count), scale = params[2])


mdl_rm2 <- lm(
  ViolentCrimesPerPop ~ .,
  data[-c(463, 1583, 1271, 669), which(colSums(is.na(data)) == 0)] |>
    dplyr::select(-communityname, -state, -fold)
)

coef_orig <- summary(mdl)$coefficients["racePctBlack",1]
coef_part1 <- summary(mdl_rm)$coefficients["racePctBlack",1]
coef_all <- summary(mdl_rm2)$coefficients["racePctBlack",1]

round(coef_all - coef_orig, 4)
round(coef_part1 - coef_orig, 4)
round(coef_all - coef_part1, 4)

coef_all/coef_orig
coef_part1/coef_orig
1/(coef_all/coef_part1)
