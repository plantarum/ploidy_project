# Front Matter - PACKAGES
#raster for projecting coordinates
library(raster)
#dismo for maxent models
library(dismo)
#for B1 and B2
library(ENMTools)

setwd("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project")
#READING DATA INTO R
wc <- getData(name = "worldclim", var = "bio", res = 5)

#parental species
parent1 <- read.csv("Datasets/Aegilops_comosa.csv")
parent1 <- parent1[parent1$decimalLongitude != 0 & parent1$decimalLatitude != 0,]
parent1 <- parent1[complete.cases(parent1$decimalLatitude),]
parent1C <- parent1[, c("decimalLongitude", "decimalLatitude")]
parent1C <- as.data.frame(subset(parent1C, !is.na(decimalLongitude)))

parent1C <- SpatialPoints(parent1C, proj4string=CRS(proj4string(wrld_simpl)))
parent1C <- parent1C[wrld_simpl] #remove points not on land
parent1C <- data.frame(coordinates(parent1C))

parent2 <- read.csv("Datasets/Aegilops_umbellulata.csv")
parent2 <- parent2[parent2$decimalLongitude != 0 & parent2$decimalLatitude != 0, ] 
parent2 <- parent2[complete.cases(parent2$decimalLatitude),]
parent2C <- parent2[, c("decimalLongitude", "decimalLatitude")]
parent2C <- as.data.frame(subset(parent2C, !is.na(decimalLongitude)))

parent2C <- SpatialPoints(parent2C, proj4string=CRS(proj4string(wrld_simpl)))
parent2C <- parent2C[wrld_simpl] #remove points not on land
parent2C <- data.frame(coordinates(parent2C))

#allopolyploid 
hybrid <- read.csv("Datasets/Aegilops_kotschyi.csv")
hybrid <- hybrid[hybrid$decimalLongitude != 0 & hybrid$decimalLatitude != 0, ]  
hybrid <- hybrid[complete.cases(hybrid$decimalLatitude),]
hybridC <- hybrid[, c("decimalLongitude", "decimalLatitude")] #only x and y coordinates 
hybridC <- as.data.frame(subset(hybridC, !is.na(decimalLongitude)))

hybridC <- SpatialPoints(hybridC, proj4string=CRS(proj4string(wrld_simpl)))
hybridC <- hybridC[wrld_simpl] #remove points not on land
hybridC <- data.frame(coordinates(hybridC)) #extract coordinates and convert to dataframe

#get name of species for saving maxent model
parent1_name <- parent1$species[1]
parent1_name <- unlist(strsplit(parent1_name, " "))
parent2_name <- parent2$species[1]
parent2_name <- unlist(strsplit(parent2_name, " "))
hybrid_name <- hybrid$species[1]
hybrid_name <- unlist(strsplit(hybrid_name, " "))

#SAMPLING BIAS 
#accounting for sampling bias using grid method

sample.bias <- function(coords){
  r <- raster(extent(range(coords[,1]), range(coords[,2]))) #create an empty raster
  res(r) <- 1/60 
  r <- extend(r, extent(r) + 1) #extend it outwards one cell
  coords <- gridSample(as.data.frame(coords), r, n = 1)
}

hybridsel <- sample.bias(hybridC) 
parent1sel <- sample.bias(parent1C)
parent2sel <- sample.bias(parent2C)


#MaxEnt Modeling 
##hybrid
xm <- maxent(wc, hybridsel) 
pr <- predict(wc, xm) 
plot(pr) #model of the species 
plot(xm) #look at how each variable responds to the species distribution

raster.breadth(pr) #get B1 and B2 suitability values 

save(xm, pr, file = paste("Models/", hybrid_name[1], "_", hybrid_name[2], "_19var", ".RData", sep= ""))
load(paste("Models/", hybrid_name[1], "_", hybrid_name[2], "_19var", ".RData", sep= ""))

##parent 1
xm_parent1 <- maxent(wc, parent1sel) 
pr_parent1 <- predict(wc, xm_parent1) 
plot(pr_parent1) #model of the species 
plot(xm_parent1) #look at how each variable responds to the species distribution

raster.breadth(pr_parent1)

save(xm_parent1, pr_parent1, file = paste("Models/", parent1_name[1], "_", parent1_name[2], "_19var",
                                          ".RData", sep= ""))
load(paste("Models/", parent1_name[1], "_", parent1_name[2], "_19var",".RData", sep= ""))


##parent2
xm_parent2 <- maxent(wc, parent2sel) 
pr_parent2 <- predict(wc, xm_parent2) 
plot(pr_parent2) #model of the species 
plot(xm_parent2) #look at how each variable responds to the species distribution

raster.breadth(pr_parent2)

save(xm_parent2, pr_parent2, file = paste("Models/", parent2_name[1], "_", parent2_name[2], "_19var",
                                           ".RData", sep= ""))
load(paste("Models/", parent2_name[1], "_", parent2_name[2], "_19var",".RData", sep= ""))

