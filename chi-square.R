## Chi-square analysis ##

setwd("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project")
#data <- read.csv("Ploidy_project_analysis_Mon_Apr_04_2022.csv")
data <- read.csv("Results_summary_onlydiploids_April2022.csv")
#model_data <- read.csv("maxentmodels_B1+B2_April2022.csv")
model_data <- read.csv("maxentmodels_B1+B2_onlydiploidPs_April2022.csv")

## Geographic range area:
ranges <- data[, c("Range_km2", "P1_range_km2", "P2_range_km2")]

## Env range area:
env_ranges <- data[, c("Env_range_area", "Env_range_area_P1", "Env_range_area_P2")]

## convert to ranks, 1 is the smallest range
rangeRanks <- apply(ranges, 1, function(x) rank(x, na.last = "keep"))
rangeRanks <- t(rangeRanks)

## Two parents known:
rangeRanks2P <- subset(rangeRanks, !is.na(rangeRanks[, "P2_range_km2"]))

## One parent known:
rangeRanks1P <- subset(rangeRanks, is.na(rangeRanks[, "P2_range_km2"]))

## data:
table(rangeRanks1P[, "Range_km2"])
table(rangeRanks2P[, "Range_km2"])

## Tests:
chisq.test(table(rangeRanks1P[, "Range_km2"]))
chisq.test(table(rangeRanks2P[, "Range_km2"]))


## Environmental Range :

envrangeRanks <- apply(env_ranges, 1, function(x) rank(x, na.last = "keep"))
envrangeRanks <- t(envrangeRanks)
envrangeRanks2P <- subset(envrangeRanks, !is.na(envrangeRanks[, "Env_range_area_P2"]))
envrangeRanks1P <- subset(envrangeRanks, is.na(envrangeRanks[, "Env_range_area_P2"]))

## data:
table(envrangeRanks1P[, "Env_range_area"])
table(envrangeRanks2P[, "Env_range_area"])

## Tests:
chisq.test(table(envrangeRanks1P[, "Env_range_area"]))
chisq.test(table(envrangeRanks2P[, "Env_range_area"]))

## Max. Latitude:
max_lat <- data[, c("Max_Lat_hybrid", "Max_Lat_P1", "Max_Lat_P2")]

maxlatRanks <- apply(max_lat, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
maxlatRanks <- t(maxlatRanks)
maxlatRanks2P <- subset(maxlatRanks, !is.na(maxlatRanks[, "Max_Lat_P2"]))
maxlatRanks1P <- subset(maxlatRanks, is.na(maxlatRanks[, "Max_Lat_P2"]))

## data:
table(maxlatRanks1P[, "Max_Lat_hybrid"])
table(maxlatRanks2P[, "Max_Lat_hybrid"])

## Tests:
chisq.test(table(maxlatRanks1P[, "Max_Lat_hybrid"]))
chisq.test(table(maxlatRanks2P[, "Max_Lat_hybrid"]))

## Max. Temperature:
max_temp <- data[, c("Temp_max", "P1_Temp_max", "P2_Temp_max")]

maxtempRanks <- apply(max_temp, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
maxtempRanks <- t(maxtempRanks)
maxtempRanks2P <- subset(maxtempRanks, !is.na(maxtempRanks[, "P2_Temp_max"]))
maxtempRanks1P <- subset(maxtempRanks, is.na(maxtempRanks[, "P2_Temp_max"]))

## data:
table(maxtempRanks1P[, "Temp_max"])
table(maxtempRanks2P[, "Temp_max"])

## Tests:
chisq.test(table(maxtempRanks1P[, "Temp_max"]))
chisq.test(table(maxtempRanks2P[, "Temp_max"]))

## Min. Temperature:
min_temp <- data[, c("Temp_min", "P1_temp_min", "P2_temp_min")]

mintempRanks <- apply(min_temp, 1, function(x) rank(x, na.last = "keep",ties.method = "last"))
mintempRanks <- t(mintempRanks)
mintempRanks2P <- subset(mintempRanks, !is.na(mintempRanks[, "P2_temp_min"]))
mintempRanks1P <- subset(mintempRanks, is.na(mintempRanks[, "P2_temp_min"]))

## data:
table(mintempRanks1P[, "Temp_min"])
table(mintempRanks2P[, "Temp_min"])

## Tests:
chisq.test(table(mintempRanks1P[, "Temp_min"]))
chisq.test(table(mintempRanks2P[, "Temp_min"]))

## Temperature breadth:
tempbreadth <- data[, c("Temp_breadth", "P1_temp_breadth", "P2_temp_breadth")]

tempbreadthRanks <- apply(tempbreadth, 1, function(x) rank(x, na.last = "keep"))
tempbreadthRanks <- t(tempbreadthRanks)
tempbreadthRanks2P <- subset(tempbreadthRanks, !is.na(tempbreadthRanks[, "P2_temp_breadth"]))
tempbreadthRanks1P <- subset(tempbreadthRanks, is.na(tempbreadthRanks[, "P2_temp_breadth"]))

## data:
table(tempbreadthRanks1P[, "Temp_breadth"])
table(tempbreadthRanks2P[, "Temp_breadth"])

## Tests:
chisq.test(table(tempbreadthRanks1P[, "Temp_breadth"]))
chisq.test(table(tempbreadthRanks2P[, "Temp_breadth"]))

## Max. Precipitation:
max_precip <- data[, c("Precip_max", "P1_precip_max", "P2_precip_max")]

maxprecipRanks <- apply(max_precip, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
maxprecipRanks <- t(maxprecipRanks)
maxprecipRanks2P <- subset(maxprecipRanks, !is.na(maxprecipRanks[, "P2_precip_max"]))
maxprecipRanks1P <- subset(maxprecipRanks, is.na(maxprecipRanks[, "P2_precip_max"]))

## data:
table(maxprecipRanks1P[, "Precip_max"])
table(maxprecipRanks2P[, "Precip_max"])

## Tests:
chisq.test(table(maxprecipRanks1P[, "Precip_max"]))
chisq.test(table(maxprecipRanks2P[, "Precip_max"]))

## Min. Precipitation:
min_precip <- data[, c("Precip_min", "P1_precip_min", "P2_precip_min")]

minprecipRanks <- apply(min_precip, 1, function(x) rank(x, na.last = "keep", ties.method = "last"))
minprecipRanks <- t(minprecipRanks)
minprecipRanks2P <- subset(minprecipRanks, !is.na(minprecipRanks[, "P2_precip_min"]))
minprecipRanks1P <- subset(minprecipRanks, is.na(minprecipRanks[, "P2_precip_min"]))

## data:
table(minprecipRanks1P[, "Precip_min"])
table(minprecipRanks2P[, "Precip_min"])

## Tests:
chisq.test(table(minprecipRanks1P[, "Precip_min"]))
chisq.test(table(minprecipRanks2P[, "Precip_min"]))

## Precipitation breadth:
precipbreadth <- data[, c("Precip_breadth", "P1_Precip_breadth", "P2_precip_breadth")]

precipbreadthRanks <- apply(precipbreadth, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
precipbreadthRanks <- t(precipbreadthRanks)
precipbreadthRanks2P <- subset(precipbreadthRanks, !is.na(precipbreadthRanks[, "P2_precip_breadth"]))
precipbreadthRanks1P <- subset(precipbreadthRanks, is.na(precipbreadthRanks[, "P2_precip_breadth"]))

## data:
table(precipbreadthRanks1P[, "Precip_breadth"])
table(precipbreadthRanks2P[, "Precip_breadth"])

## Tests:
chisq.test(table(precipbreadthRanks1P[, "Precip_breadth"]))
chisq.test(table(precipbreadthRanks2P[, "Precip_breadth"]))

## Parents vs Parent-hybrid overlap for environmental range:
overlap <- data[, c("P1_hybrid_overlap", "P2_hybrid_overlap", "Parents_env_overlap")]

overlapRanks <- apply(overlap, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
overlapRanks <- t(overlapRanks)
overlapRanks2P <- subset(overlapRanks, !is.na(overlapRanks[, "P2_hybrid_overlap"]))

## data:
table(overlapRanks2P[, "Parents_env_overlap"])

## Tests:
chisq.test(table(overlapRanks2P[, "Parents_env_overlap"]))

## Parents vs Parent-hybrid overlap for geographical range:
geo_overlap <- data[, c("P1_h_geo_overlap", "P2_h_geo_overlap", "Parents_geo_overlap")]

geo_overlapRanks <- apply(geo_overlap, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
geo_overlapRanks <- t(geo_overlapRanks)
geo_overlapRanks2P <- subset(geo_overlapRanks, !is.na(geo_overlapRanks[, "P2_h_geo_overlap"]))

## data:
table(geo_overlapRanks2P[, "Parents_geo_overlap"])

## Tests:
chisq.test(table(geo_overlapRanks2P[, "Parents_geo_overlap"]))

## Parents vs Parent-hybrid overlap for B1+B2:
B1 <- model_data[, c("B1", "B1_P1", "B1_P2")]

B1Ranks <- apply(B1, 1, function(x) rank(x, na.last = "keep"))
B1Ranks <- t(B1Ranks)
B1Ranks2P <- subset(B1Ranks, !is.na(B1Ranks[, "B1_P2"]))
B1Ranks1P <- subset(B1Ranks, is.na(B1Ranks[, "B1_P2"]))

## data:
table(B1Ranks1P[, "B1"])
table(B1Ranks2P[, "B1"])

## Tests:
chisq.test(table(B1Ranks1P[, "B1"]))
chisq.test(table(B1Ranks2P[, "B1"]))

B2 <- model_data[, c("B2", "B2_P1", "B2_P2")]

B2Ranks <- apply(B2, 1, function(x) rank(x, na.last = "keep"))
B2Ranks <- t(B2Ranks)
B2Ranks2P <- subset(B2Ranks, !is.na(B2Ranks[, "B2_P2"]))
B2Ranks1P <- subset(B2Ranks, is.na(B2Ranks[, "B2_P2"]))

## data:
table(B2Ranks1P[, "B2"])
table(B2Ranks2P[, "B2"])

## Tests:
chisq.test(table(B2Ranks1P[, "B2"]))
chisq.test(table(B2Ranks2P[, "B2"]))

## Centroid position:
centroid <- data[, c("ABS.centroid_latitude.", "ABS.P1_centroid_latitude.", "ABS.P2_centroid_latitude.")]

centroidRanks <- apply(centroid, 1, function(x) rank(x, na.last = "keep", ties.method = "first"))
centroidRanks <- t(centroidRanks)
centroidRanks2P <- subset(centroidRanks, !is.na(centroidRanks[, "ABS.P2_centroid_latitude."]))
centroidRanks1P <- subset(centroidRanks, is.na(centroidRanks[, "ABS.P2_centroid_latitude."]))

## data:
table(centroidRanks1P[, "ABS.centroid_latitude."])
table(centroidRanks2P[, "ABS.centroid_latitude."])

## Tests:
chisq.test(table(centroidRanks1P[, "ABS.centroid_latitude."]))
chisq.test(table(centroidRanks2P[, "ABS.centroid_latitude."]))
