setwd("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project")

data <- read.csv("result_summary.csv")
sdm_data <- read.csv("maxentmodels_B1+B2.csv")

#data without non diploid parents 
summ <- read.csv("maxentmodels_B1+B2_onlydiploidPs.csv")
data2 <- read.csv("result_summary_onlydiploidPs.csv")

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#RANDOMIZATION TESTS
#necessary functions
permuteRanks <- function(x){
  ## Permutes the ranks of a vector, ignoring NA values
  ## Should be able to correctly handle vectors of any length
  maxSize <- length(x)
  Values <- x[!is.na(x)]
  groupSize <- sum(!is.na(x))
  
  perm <- sample(Values)
  if(groupSize < maxSize)
    perm <- c(perm, rep(NA, maxSize - groupSize))
  return(perm)
}  

rankDist <- function(x, n = 1000){
  res = numeric(n)
  for(i in 1:n){
    permutedRanks <- t(apply(x, 1, permuteRanks))
    res[i] <- mean(permutedRanks[ , 1])
  }
  return(res)
}

#------------------------------------------------------
#Ranges randomization tests
## for convenience pull out the range columns
ranges <- data[, c("Range_km2", "P1_range_km2", "P2_range_km2", "P3_range_km2")]
env_ranges <- data[, c("Env_range", "P1_Env_range", "P2_Env_range", "P3_Env_range" )]

## convert to ranks, 1 is the smallest range
rangeRanks <- apply(ranges, 1, function(x) rank(x, na.last = "keep"))
rangeRanks <- t(rangeRanks)

envrangeRanks <- apply(env_ranges, 1, function(x) rank(x, na.last = "keep"))
envrangeRanks <- t(envrangeRanks)

## Calculate the NULL distribution for geographical range:

res <- rankDist(rangeRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid geographic rank",
     sub = "Larger Rank = Larger Range",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(rangeRanks[, "Range_km2"]))

#p-value
sum(res < mean(rangeRanks[, "Range_km2"]))/length(res) * 2

## Calculate the NULL distribution for environmental range:

res <- rankDist(envrangeRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid environment rank",  sub = "Larger Rank = Larger Range",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(envrangeRanks[, "Env_range"]))

#p-value
sum(res < mean(envrangeRanks[, "Env_range"]))/length(res) * 2

#---------------------------------------------------------------
#Range Randomization tests without parents that aren't diploid
ranges <- data2[, c("Range_km2", "P1_range_km2", "P2_range_km2")]
env_ranges <- data2[, c("Env_range", "P1_Env_range", "P2_Env_range")]

## convert to ranks, 1 is the smallest range
rangeRanks <- apply(ranges, 1, function(x) rank(x, na.last = "keep"))
rangeRanks <- t(rangeRanks)

envrangeRanks <- apply(env_ranges, 1, function(x) rank(x, na.last = "keep"))
envrangeRanks <- t(envrangeRanks)

## Calculate the NULL distribution for geographical range:

res <- rankDist(rangeRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid geographic rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(rangeRanks[, "Range_km2"]))

#p-value
sum(res < mean(rangeRanks[, "Range_km2"]))/length(res) * 2

## Calculate the NULL distribution for environmental range:

res <- rankDist(envrangeRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid environment rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(envrangeRanks[, "Env_range"]))

#p-value
sum(res < mean(envrangeRanks[, "Env_range"]))/length(res) * 2
#---------------------------------------------------------------
##B1 and B2 randomization tests

B1s <- sdm_data[, c("B1", "B1_P1", "B1_P2", "B1_P3")]
B2s <- sdm_data[, c("B2", "B2_P1", "B2_P2", "B2_P3")]

## convert to ranks, 1 is the smallest 
B1sRanks <- apply(B1s, 1, function(x) rank(x, na.last = "keep"))
B1sRanks <- t(B1sRanks)

B2sRanks <- apply(B2s, 1, function(x) rank(x, na.last = "keep"))
B2sRanks <- t(B2sRanks)

## Calculate the NULL distribution for B1:
res <- rankDist(B1sRanks, n = 10000)

#p-value
sum(res < mean(B1sRanks[, "B1"]))/length(res) * 2

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid B1 rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(B1sRanks[, "B1"]))

##Calculate the NULL distribution for B2:
res <- rankDist(B2sRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid B2 rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(B2sRanks[, "B2"]))

#p-value
sum(res < mean(B2sRanks[, "B2"]))/length(res) * 2

#----------------------------------------------------------------------
##B1 and B2 with only parents that are diploid
##note: do not need to include P3 since one of the parents isn't diploid for the species with 3 species 
B1s <- summ[, c("B1", "B1_P1", "B1_P2")] 
B2s <- summ[, c("B2", "B2_P1", "B2_P2")]

## convert to ranks, 1 is the smallest 
B1sRanks <- apply(B1s, 1, function(x) rank(x, na.last = "keep"))
B1sRanks <- t(B1sRanks)

B2sRanks <- apply(B2s, 1, function(x) rank(x, na.last = "keep"))
B2sRanks <- t(B2sRanks)

## Calculate the NULL distribution for B1:

res <- rankDist(B1sRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid B1 rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(B1sRanks[, "B1"]))

#p-value
sum(res < mean(B1sRanks[, "B1"]))/length(res) * 2

##Calculate the NULL distribution for B2:
res <- rankDist(B2sRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid B2 rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(B2sRanks[, "B2"]))

#p-value
sum(res < mean(B2sRanks[, "B2"]))/length(res) * 2

#------------------------------------------------------------------------
##Climate variables (temp_max, temp_min, temp_breadth, precip_max, precip_min, precip_seasonality)

temp_max <- data[, c("Temp_max","P1_Temp_max","P2_Temp_max", "P3_Temp_max")]
temp_min <- data[, c("Temp_min","P1_temp_min","P2_temp_min", "P3_temp_min")]
temp_breadth <- data[, c("Temp_breadth","P1_temp_breadth","P2_temp_breadth", "P3_temp_breadth")]
precip_max <- data[, c("Precip_max","P1_precip_max","P2_precip_max", "P3_precip_max")]
precip_min <- data[, c("Precip_min","P1_precip_min","P2_precip_min", "P3_precip_min")]
precip_seasonality <- data[, c("Precip_seasonality","P1_precip_seasonality","P2_precip_seasonality",
                               "P3_precip_seasonality")]

## convert to ranks, 1 is the smallest 
tempmaxRanks <- apply(temp_max, 1, function(x) rank(x, na.last = "keep"))
tempmaxRanks <- t(tempmaxRanks)

tempminRanks <- apply(temp_min, 1, function(x) rank(x, na.last = "keep"))
tempminRanks <- t(tempminRanks)

tempbreadthRanks <- apply(temp_breadth, 1, function(x) rank(x, na.last = "keep"))
tempbreadthRanks <- t(tempbreadthRanks)

precipmaxRanks <- apply(precip_max, 1, function(x) rank(x, na.last = "keep"))
precipmaxRanks <- t(precipmaxRanks)

precipminRanks <- apply(precip_min, 1, function(x) rank(x, na.last = "keep"))
precipminRanks <- t(precipminRanks)

precipSRanks <- apply(precip_seasonality, 1, function(x) rank(x, na.last = "keep"))
precipSRanks <- t(precipSRanks)

## Calculate the NULL distribution for temp_max:

res <- rankDist(tempmaxRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp max rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempmaxRanks[, "Temp_max"]))

#p-value
sum(res > mean(tempmaxRanks[, "Temp_max"]))/length(res) * 2

##Calculate the NULL distribution for temp_min
res <- rankDist(tempminRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp min rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempminRanks[, "Temp_min"]))

#p-value
sum(res > mean(tempminRanks[, "Temp_min"]))/length(res) * 2

##Calculate the NULL distribution for temp_breadth:
res <- rankDist(tempbreadthRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp breadth rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempbreadthRanks[, "Temp_breadth"]))

#p-value
sum(res < mean(tempbreadthRanks[, "Temp_breadth"]))/length(res) * 2

## Calculate the NULL distribution for precip_max:

res <- rankDist(precipmaxRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip max rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipmaxRanks[, "Precip_max"]))

#p-value
sum(res < mean(precipmaxRanks[, "Precip_max"]))/length(res) * 2

##Calculate the NULL distribution for precip_min:
res <- rankDist(precipminRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip min rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipminRanks[, "Precip_min"]))

#p-value
sum(res > mean(precipminRanks[, "Precip_min"]))/length(res) * 2

##Calculate the NULL distribution for precip_seasonality:
res <- rankDist(precipSRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip seasonality rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipSRanks[, "Precip_seasonality"]))

#p-value
sum(res < mean(precipSRanks[, "Precip_seasonality"]))/length(res) * 2

#------------------------------------------------------------------------
##Climate variables (temp_max, temp_min, temp_breadth, precip_max, precip_min, precip_seasonality) w/o 
#parents that aren't diploid 

temp_max <- data2[, c("Temp_max","P1_Temp_max","P2_Temp_max")]
temp_min <- data2[, c("Temp_min","P1_temp_min","P2_temp_min")]
temp_breadth <- data2[, c("Temp_breadth","P1_temp_breadth","P2_temp_breadth")]
precip_max <- data2[, c("Precip_max","P1_precip_max","P2_precip_max")]
precip_min <- data2[, c("Precip_min","P1_precip_min","P2_precip_min")]
precip_seasonality <- data2[, c("Precip_seasonality","P1_precip_seasonality","P2_precip_seasonality")]

## convert to ranks, 1 is the smallest 
tempmaxRanks <- apply(temp_max, 1, function(x) rank(x, na.last = "keep"))
tempmaxRanks <- t(tempmaxRanks)

tempminRanks <- apply(temp_min, 1, function(x) rank(x, na.last = "keep"))
tempminRanks <- t(tempminRanks)

tempbreadthRanks <- apply(temp_breadth, 1, function(x) rank(x, na.last = "keep"))
tempbreadthRanks <- t(tempbreadthRanks)

precipmaxRanks <- apply(precip_max, 1, function(x) rank(x, na.last = "keep"))
precipmaxRanks <- t(precipmaxRanks)

precipminRanks <- apply(precip_min, 1, function(x) rank(x, na.last = "keep"))
precipminRanks <- t(precipminRanks)

precipSRanks <- apply(precip_seasonality, 1, function(x) rank(x, na.last = "keep"))
precipSRanks <- t(precipSRanks)

## Calculate the NULL distribution for temp_max:

res <- rankDist(tempmaxRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp max rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempmaxRanks[, "Temp_max"]))

#p-value
sum(res > mean(tempmaxRanks[, "Temp_max"]))/length(res) * 2

##Calculate the NULL distribution for temp_min
res <- rankDist(tempminRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp min rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempminRanks[, "Temp_min"]))

#p-value
sum(res > mean(tempminRanks[, "Temp_min"]))/length(res) * 2

##Calculate the NULL distribution for temp_breadth:
res <- rankDist(tempbreadthRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid temp breadth rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(tempbreadthRanks[, "Temp_breadth"]))

#p-value
sum(res < mean(tempbreadthRanks[, "Temp_breadth"]))/length(res) * 2

## Calculate the NULL distribution for precip_max:

res <- rankDist(precipmaxRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip max rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipmaxRanks[, "Precip_max"]))

#p-value
sum(res < mean(precipmaxRanks[, "Precip_max"]))/length(res) * 2

##Calculate the NULL distribution for precip_min:
res <- rankDist(precipminRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip min rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipminRanks[, "Precip_min"]))

#p-value
sum(res > mean(precipminRanks[, "Precip_min"]))/length(res) * 2

##Calculate the NULL distribution for precip_seasonality:
res <- rankDist(precipSRanks, n = 10000)

hist(res, breaks = 20,
     main = "Null distribution of mean polyploid precip seasonality rank",
     xlab = "Polyploid Rank")
## 95% confidence interval, in red:
abline(v = quantile(res, probs = c(0.025, 0.975)), col = 2)
## Actual mean rank, in black:
abline(v = mean(precipSRanks[, "Precip_seasonality"]))

#p-value
sum(res < mean(precipSRanks[, "Precip_seasonality"]))/length(res) * 2
