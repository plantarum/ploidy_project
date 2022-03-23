#Clear Environment
rm(list=ls())

############################################ PROGRAM DESCRIPTION
# The purpose of this program is to analyze the occurrence data from GBIF for allopolyploid and their parental species

############################################ OUTPUT FILES
# A .csv file with data from analysis following this trend : "Ploidy_project_analysis", date.str , ".csv" 

#-------------------------------------------------------------------------------------------------------------------------------------
# Front Matter - PACKAGES
# maptools for mapping ranges
library(maptools)
library(raster)
library(ENMTools)
library(dismo)
library(ade4)

#-------------------------------------------------------------------------------------------------------------------------------------
# Session Tool Information: 
options(width = 100)
devtools::session_info()

#-------------------------------------------------------------------------------------------------------------------------------------
# Front Matter - PATHS

#-------------------------------------------------------------------------------------------------------------------------------------
# Front Matter - FUNCTIONS
# test.MBplots function used to test each file for white noise using normwhn.test's whitenoise.test function and Box.test

#-------------------------------------------------------------------------------------------------------------------------------------
# Front Matter - DATA
# In path: 
#Occurrence data downloaded from GBIF: e.g. Aegilops_cylindrica.txt
#Worldclim data from www.worldclim.org
#wrld_simpl map

#-------------------------------------------------------------------------------------------------------------------------------------
######################################################################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------
# Begin Program Area
#-------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#setwd("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project")

#read in information about the hybrid and its parents 
comb <- read.csv("Hybrid_parent_comb.csv") #csv with hybrid and its parents
hybrids <- comb[,1]
hybrids <- gsub(" ", "", hybrids) #remove extra spaces
parent1s <- comb[,2]
parent1s <- gsub(" ", "", parent1s)
parent2s <- comb[,3]
parent2s <- gsub(" ", "", parent2s)

#set date string
if (length(unlist(strsplit(date(), " "))) == 5) date.str <- paste(unlist(strsplit(date(), " "))[c(1:3,5)], collapse='_')
if (length(unlist(strsplit(date(), " "))) == 6) date.str <- paste(unlist(strsplit(date(), " "))[c(1:2,4,6)], collapse='_') 

#get worldclim data 
wc <- getData(name = "worldclim", var = "bio", res = 10, path = "../")
res_km <- res(wc) *111

#convert worldclim data to 10000 randomly sampled points across the Earth
Env <- rasterToPoints(wc$bio1)
Env <- Env[,1:2] #extract xy points for Earth 

glob_geo <- Env[sample(nrow(Env), 10000), ]

#function to get the coordinates from the csv files
getcoords <- function(csv){
  coords <- csv[csv$decimalLongitude !=0 & csv$decimalLatitude !=0, ] #remove zero coordinates
  coords <- subset(coords, !is.na(decimalLongitude)) #remove NAs
  coords <- coords[, c("decimalLongitude", "decimalLatitude")] #just the coordinates
}

#function to extract the climatic conditions where the species occur
envextract <- function(sp.coords){
  env <- extract(wc, sp.coords)
  env <- env[complete.cases(env), ]
}


#function to calculate min, max and temperature breadth for where the species occurs
tempbreadth <- function(sp.coords){
  bio5 <- extract(wc$bio5, sp.coords)
  maxtemp <- max(bio5, na.rm=T)
  
  bio6 <- extract(wc$bio6, sp.coords)
  mintemp <- min(bio6, na.rm=T)
  
  tempbreadth <- maxtemp - mintemp
  
  results <- c(maxtemp, mintemp, tempbreadth)
  return(results)
}

precipbreadth <- function(sp.coords){
  Bio13<- extract(wc$bio13, sp.coords) 
  maxprec <- max(Bio13, na.rm=T)
  
  Bio14 <- extract(wc$bio14, sp.coords)
  minprec <- min(Bio14, na.rm=T)
  
  precbreadth <- extract(wc$bio15, sp.coords)
  precbreadth <- max(precbreadth, na.rm = T)
  
  results <- c(maxprec, minprec, precbreadth)
}

env.range.area <- function(pca.scores){
  mask.rast <- raster(extent(range(pca.scores[,1]), range(pca.scores[,2]))) #create an empty raster
  mask.rast <- extend(mask.rast, extent(mask.rast) + 1)
  res(mask.rast) <- 0.2
  
  #set the background cells in the raster to 0
  mask.rast[!is.na(mask.rast)] <- 0 
  
  #set the cells that contain points to 1 
  env_range <- rasterize(pca.scores, mask.rast, field = 1)
  env_range <- merge(env_range, mask.rast)
  
  #determine resolution of the raster
  res_rast <- res(env_range)
  
  #determine number of cells that have value of 1
  ncells <- freq(env_range, value= 1, useNA= "no")
  
  #calculate area
  area <- res_rast[1] * res_rast[2] * ncells
}

geo.range.area <- function(coordinates){
  #Calculating Area range and overlap for the species
  ##create enmtool species objects
  enm <- enmtools.species(presence.points = coordinates)
  
  ##create a raster layer with points with 100km buffer
  enm$range <- background.raster.buffer(enm$presence.points, 100000, wc)
  
  ##Calculate the area of the range
  #ncells <- freq(hybrid_enm$range, value = 1, useNA= "no") #number of cells where the hybrid is present
  ncells <- freq(enm$range, value = 1, useNA= "no") #number of cells where the hybrid is present
  area <- res_km[1] * res_km[2] * ncells #range size is equal to the resolution times the number of cells
  #area <- res_km[1] *res_km[2] * ncells_hybrid #range size is equal to the resolution times the number of cells
}



data("wrld_simpl")

Summary.File <- NULL
i <- 1

Aegilops_cylindrica.csv; Aegilops_tauschii.csv; Aegilops_caudata.csv

for (i in 1:nrow(comb)){
  #get the names of the species
  hybrid_name <- hybrids[i]
  parent1_name <- parent1s[i]

  hybrid_name <- "Aegilops_cylindrica"
  parent1_name <- "Aegilops_tauschii"
  parent2_name <- "Aegilops_caudata"

  #read in data for the species and extract coordinates 
  #hybrid <- read.csv(paste("../Datasets/", hybrid_name, ".csv", sep=""))
  hybrid <- read.csv(paste("Datasets/", hybrids[i], ".csv", sep="")) 
  hybridC <- getcoords(hybrid)
  
  #parent1 <- read.csv(paste("../Datasets/", parent1_name, ".csv", sep=""))
  parent1 <- read.csv(paste("Datasets/", parent1s[i], ".csv", sep=""))
  parent1C <- getcoords(parent1)

  #determining max,min temp and precipitation breadths for each species and parent 1 
  hybrid_temptest <- tempbreadth(hybridC)
  hybrid_tempmax <- hybrid_temptest[1]
  hybrid_tempmin <-hybrid_temptest[2]
  hybrid_tempbreadth <- hybrid_temptest[3]
  
  parent1_temptest <- tempbreadth(parent1C)
  parent1_tempmax <- parent1_temptest[1]
  parent1_tempmin <- parent1_temptest[2]
  parent1_tempbreadth <- parent1_temptest[3]
  
  hybrid_prectest <- precipbreadth(hybridC)
  hybrid_precmax <- hybrid_prectest[1]
  hybrid_precmin <-hybrid_prectest[2]
  hybrid_precbreadth <- hybrid_prectest[3]
  
  parent1_prectest <- precipbreadth(parent1C)
  parent1_precmax <- parent1_prectest[1]
  parent1_precmin <- parent1_prectest[2]
  parent1_precbreadth <- parent1_prectest[3]

  hybrid_enm <- enmtools.species(presence.points = hybridC)
  hybrid_enm$range <- background.raster.buffer(hybrid_enm$presence.points,
                                               100000, wc) 
  parent1_enm <- enmtools.species(presence.points = parent1C)
  parent1_enm$range <- background.raster.buffer(parent1_enm$presence.points,
                                               100000, wc) 

  #Geographical ranges
  area_hybrid <- geo.range.area(hybridC)
  area_parent1 <- geo.range.area(parent1C)
  
  #Geographical range overlap
  comp1_overlap_geo <- geog.range.overlap(hybrid_enm, parent1_enm)
  
  #Extract climate conditions where the species occurs 
  parent1_env <- envextract(parent1C)
  
  hybrid_env <- envextract(hybridC)
  
  glob <- sampleRandom(wc, 10000) #sample a random 10000 points of environmental conditions to represent globe 
  pca.env <- dudi.pca(glob, scannf=F, nf = 2)
  ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)
  
  ### PCA scores for the whole study area
  scores.globclim <- pca.env$li
  ## PCA scores for parent1
  scores.parent1 <- suprow(pca.env, parent1_env)$li
  ## PCA scores for hybrid
  scores.hybrid <- suprow(pca.env, hybrid_env)$li
  
  #Environmental ranges
  env_range_area_hybrid<- env.range.area(scores.hybrid)
  env_range_area_P1 <- env.range.area(scores.parent1)
  
  #Environmental overlap
  grid.clim.parent1 <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                          glob1 = scores.globclim,
                                          sp = scores.parent1, R=100,
                                          th.sp=0)
  
  
  grid.clim.hybrid <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                           glob1 = scores.globclim,
                                           sp = scores.hybrid, R=100,
                                           th.sp=0)
  
  D.overlap <- ecospat.niche.overlap(grid.clim.parent1, grid.clim.hybrid,
                                     cor = TRUE)
  Schoeners_metric <- D.overlap$D
  Hellinger_metric <- D.overlap$I
  
  jpeg(filename = paste("Plots/niche_overlap_for_", hybrid_name, "+_", parent1_name,".jpeg", sep= ""))
  
  ecospat.plot.niche.dyn(grid.clim.parent1, grid.clim.hybrid, quant = 0.25,
                         interest = 2, title = "Niche Overlap",
                         name.axis1 = "PC1", name.axis2 = "PC2")
  dev.off() #prevents the display of the plots 
  
  index_parent1comp <- ecospat.niche.dyn.index(grid.clim.parent1, grid.clim.hybrid) #determining niche overlap b/w parent 1 and hybrid
  comp1_overlap <- index_parent1comp$dynamic.index.w["stability"]
  comp1_hybrid_only <- index_parent1comp$dynamic.index.w["expansion"]
  parent1_only <- index_parent1comp$dynamic.index.w["unfilling"]
  
  #If the hybrid has a parent 2 then redo the above calculations with parent 2:
  if (parent2s[i] != "") {
    #parent2_name <- "Aegilops_caudata"
    parent2_name <- parent2s[i]
    #parent2 <- read.csv(paste("../Datasets/", parent2_name, ".csv", sep=""))
    parent2 <- read.csv(paste("../Datasets/", parent2s[i], ".csv", sep=""))
    parent2C <- getcoords(parent2)
    
    #max, min and breadths for temperature and precipitation
    parent2_temptest <- tempbreadth(parent2C)
    parent2_tempmax <- parent2_temptest[1]
    parent2_tempmin <- parent2_temptest[2]
    parent2_tempbreadth <- parent2_temptest[3]
    
    parent2_prectest <- precipbreadth(parent2C)
    parent2_precmax <- parent2_prectest[1]
    parent2_precmin <- parent2_prectest[2]
    parent2_precbreadth <- parent2_prectest[3]

    parent2_enm <- enmtools.species(presence.points = parent2C)
    parent2_enm$range <- background.raster.buffer(parent2_enm$presence.points,
                                                  100000, wc) 
  
    #geographical range size and overlap
    area_parent2 <- geo.range.area(parent2C)
    comp2_overlap_geo <- geog.range.overlap(hybrid_enm, parent2_enm)
    
    #Extract environmental occurrence conditions
    parent2_env <- envextract(parent2C)
    
    ##PCA scores for parent2
    scores.parent2 <- suprow(pca.env, parent2_env)$li
    
    #Environmental range
    env_range_area_P2 <- env.range.area(scores.parent2)
    
    #Environmental range overlap 
    grid.clim.parent2 <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                               glob1 = scores.globclim,
                                               sp = scores.parent2, R=100,
                                               th.sp=0)
    
    D.overlap.2 <- ecospat.niche.overlap(grid.clim.parent2, grid.clim.hybrid,
                                       cor = TRUE)
    Schoeners_metric_2 <- D.overlap.2$D
    Hellinger_metric_2 <- D.overlap.2$I
    
    jpeg(filename = paste("Plots/Niche_overlap_for", hybrid_name, "+_", parent2_name,".jpeg", sep= ""))
    
    ecospat.plot.niche.dyn(grid.clim.parent2, grid.clim.hybrid, quant = 0.25,
                           interest = 2, title = "Niche Overlap",
                           name.axis1 = "PC1", name.axis2 = "PC2")
    dev.off()
    
    ##determining niche overlap b/w parent 2 and hybrid
    index_parent2comp <- ecospat.niche.dyn.index(grid.clim.parent2, grid.clim.hybrid)
    comp2_overlap <- index_parent2comp$dynamic.index.w["stability"]
    comp2_hybrid_only <- index_parent2comp$dynamic.index.w["expansion"]
    parent2_only <- index_parent2comp$dynamic.index.w["unfilling"]
    
  }
  
  #If the hybrid does not have a parent 2 then define all of the variables as "NA"
  if (parent2s[i] == ""){
    parent2_name <- "NA"
    area_parent2 <- "NA"
    env_range_area_P2 <- "NA"
    comp2_overlap <- "NA"
    parent2_only <- "NA"
    comp2_hybrid_only <- "NA"
    Schoeners_metric_2 <- "NA"
    Hellinger_metric_2 <- "NA"
    comp2_overlap_geo <- "NA"
    parent2_temptest <- "NA"
    parent2_tempmax <- "NA"
    parent2_tempmin <- "NA"
    parent2_tempbreadth <- "NA"
    parent2_prectest <- "NA"
    parent2_precmax <- "NA"
    parent2_precmin <- "NA"
    parent2_precbreadth <- "NA"
    
    
  }
  
  data.line <- c(hybrid_name, parent1_name, parent2_name,hybrid_tempmax,
                 parent1_tempmax, parent2_tempmax, hybrid_tempmin,
                 parent1_tempmin, parent2_tempmin, hybrid_tempbreadth,
                 parent1_tempbreadth, parent2_tempbreadth,hybrid_precmax,
                 parent1_precmax, parent2_precmax, hybrid_precmin,
                 parent1_precmin, parent2_precmin, hybrid_precbreadth,
                 parent1_precbreadth, parent2_precbreadth, area_hybrid,
                 area_parent1, comp1_overlap_geo, area_parent2,
                 comp2_overlap_geo, env_range_area_hybrid, env_range_area_P1,
                 Schoeners_metric, Hellinger_metric,env_range_area_P1,
                 Schoeners_metric_2, Hellinger_metric_2, comp1_overlap,
                 comp1_hybrid_only, parent1_only, comp2_overlap,
                 comp2_hybrid_only, parent2_only)
  Summary.File <- rbind(Summary.File, data.line)
}

column_names <- list("Species","P1", "P2", "Temp_max", "P1_Temp_max", "P2_Temp_max", "Temp_min", "P1_temp_min", "P2_temp_min", 
                     "Temp_breadth","P1_temp_breadth", "P2_temp_breadth","Precip_max","P1_precip_max", "P2_precip_max", "Precip_min",
                     "P1_precip_min","P2_precip_min", "Precip_breadth","P1_Precip_breadth","P2_precip_breadth", "Range_km2", 
                     "P1_range_km2","P1_hybrid_overlap_geo","P2_range_km2","P2_hybrid_overlap_geo","Env_range_area","Env_range_area_P1","Schoeners_overlap_metric_P1&H",
                     "Hellinger_metric_P1&H","Env_range_area_P2" ,"Schoeners_overlap_metric_P2&H","Hellinger_metric_P2&H","P1_hybrid_overlap", 
                     "Hybrid_only", "P1_only", "P2_hybrid_overlap","Hybrid_only", "P2_only")
header.names <- unlist(column_names)
header.names <- c(header.names)
Summary.File <- data.frame(Summary.File)
colnames(Summary.File) <- header.names

write.csv(Summary.File, row.names = F, file = (paste("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project/", "Ploidy_project_analysis_", date.str, ".csv", sep="")))
