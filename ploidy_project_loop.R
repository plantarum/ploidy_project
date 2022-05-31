#Clear Environment
rm(list=ls())

############################################ PROGRAM DESCRIPTION
# The purpose of this program is to analyze the occurrence data from GBIF for allopolyploid and their parental species

############################################ OUTPUT FILES
# A .csv file with data from analysis following this trend : "Ploidy_project_analysis", date.str , ".csv" 

#-------------------------------------------------------------------------------------------------------------------------------------
# Front Matter - PACKAGES
library(ecospat)
library(raster)
library(maptools)
library(dismo)
library(rgeos)
library(ENMTools)
library(geosphere)

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
data("wrld_simpl")

#-------------------------------------------------------------------------------------------------------------------------------------
######################################################################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------
# Begin Program Area
#-------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project")

#read in information about the hybrid and its parents 
comb <- read.csv("Hybrid_parent_comb.csv") #csv with hybrid and its parents
hybrids <- comb[,1]
hybrids <- gsub(" ", "", hybrids) #remove extra spaces
parent1s <- comb[,2]
parent1s <- gsub(" ", "", parent1s)
parent2s <- comb[,3]
parent2s <- gsub(" ", "", parent2s)
parent2s[parent2s == ""] <- NA

#set date string
if (length(unlist(strsplit(date(), " "))) == 5) 
  date.str <- paste(unlist(strsplit(date(), " "))[c(1:3,5)], collapse='_')
if (length(unlist(strsplit(date(), " "))) == 6) 
  date.str <- paste(unlist(strsplit(date(), " "))[c(1:2,4,6)], collapse='_') 

#get worldclim data 
wc <- getData(name = "worldclim", var = "bio", res = 10)
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
tempbreadth <- function(env){
  bio5 <- env[, "bio5"]
  maxtemp <- max(bio5, na.rm=T)
  
  bio6 <- env[, "bio6"]
  mintemp <- min(bio6, na.rm=T)
  
  tempbreadth <- maxtemp - mintemp
  
  results <- c(maxtemp, mintemp, tempbreadth)
  return(results)
}

precipbreadth <- function(env){
  Bio13 <- env[, "bio13"]
  maxprec <- max(Bio13, na.rm=T)
  
  Bio14 <- env[, "bio14"]
  minprec <- min(Bio14, na.rm=T)
  
  precbreadth <- env[, "bio15"]
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
  area <- res_rast[1] *res_rast[2] *ncells
}


Summary.File <- NULL
i <- 1

## glob <- sampleRandom(wc, 20000) #sample a random 20000 points of
##                                         #environmental conditions to represent
##                                         #globe  
glob <- getValues(wc)
glob <- glob[complete.cases(glob), ]

pca.env <- dudi.pca(glob, scannf=F, nf = 2)

## PCA scores for the whole study area
scores.globclim <- pca.env$li

for (i in 1:nrow(comb)){
  message("**** I = ", i, " ****")
  #get the names of the species
  hybrid_name <- hybrids[i]
  parent1_name <- parent1s[i]
  
  message(" Processing: ", hybrid_name)
  
  #read in data for the species and extract coordinates 
  hybrid <- read.csv(paste("Datasets/", hybrids[i], ".csv", sep="")) 
  parent1 <- read.csv(paste("Datasets/", parent1s[i], ".csv", sep=""))
  
  message("  extracting coordinates") 
  hybridC <- getcoords(hybrid)
  parent1C <- getcoords(parent1)
  
  hybridC <- SpatialPoints(hybridC, proj4string=CRS(proj4string(wrld_simpl)))
  hybridC <- hybridC[wrld_simpl] #remove points not on land
  hybridC <- data.frame(coordinates(hybridC)) #extract coordinates and convert to dataframe
  
  parent1C <- SpatialPoints(parent1C, proj4string=CRS(proj4string(wrld_simpl)))
  parent1C <- parent1C[wrld_simpl] #remove points not on land
  parent1C <- data.frame(coordinates(parent1C))
  
  #Determine max and min latitude 
  maxlat_hybrid <- max(hybridC$decimalLatitude)
  maxlat_parent1 <- max(parent1C$decimalLatitude)
  
  minlat_hybrid <- min(hybridC$decimalLatitude)
  minlat_parent1 <- min(parent1C$decimalLatitude)
  
  #Determine the centroid
  hybrid_centroid <- centroid(hybridC)
  centroid_latitude <- hybrid_centroid[2]
  centroid_longitude <- hybrid_centroid[1]
  parent1_centroid <- centroid(parent1C)
  P1_centroid_latitude <- parent1_centroid[2]
  P1_centroid_longitude <- parent1_centroid[1]
  
  ## Extract climate conditions where the species occurs 
  
  message("  extracting climate data")
  hybrid_env <- envextract(hybridC)
  parent1_env <- envextract(parent1C)
  
  #determining max,min temp and precipitation breadths for each species and parent 1 
  hybrid_temptest <- tempbreadth(hybrid_env)
  hybrid_tempmax <- hybrid_temptest[1]
  hybrid_tempmin <-hybrid_temptest[2]
  hybrid_tempbreadth <- hybrid_temptest[3]
  
  parent1_temptest <- tempbreadth(parent1_env)
  parent1_tempmax <- parent1_temptest[1]
  parent1_tempmin <- parent1_temptest[2]
  parent1_tempbreadth <- parent1_temptest[3]
  
  hybrid_prectest <- precipbreadth(hybrid_env)
  hybrid_precmax <- hybrid_prectest[1]
  hybrid_precmin <-hybrid_prectest[2]
  hybrid_precbreadth <- hybrid_prectest[3]
  
  parent1_prectest <- precipbreadth(parent1_env)
  parent1_precmax <- parent1_prectest[1]
  parent1_precmin <- parent1_prectest[2]
  parent1_precbreadth <- parent1_prectest[3]

  message("   geographic overlap")
  ## Geographical range overlap
  hybrid_enm <- enmtools.species(presence.points = hybridC)
  hybrid_enm$range <- background.raster.buffer(hybrid_enm$presence.points, 100000, wc)
  ncells <- freq(hybrid_enm$range, value = 1, useNA= "no")
  area_hybrid <- res_km[1] *res_km[2] * ncells 
  ncells <- NULL ## erase before we reuse, just in case?
  
  parent1_enm <- enmtools.species(presence.points = parent1C)
  parent1_enm$range <- background.raster.buffer(parent1_enm$presence.points, 100000, wc)
  ncells <- freq(parent1_enm$range, value = 1, useNA= "no")
  area_parent1 <- res_km[1] *res_km[2] * ncells 
  
  message("   geographic overlap")
  comp1_overlap_geo <- geog.range.overlap(hybrid_enm, parent1_enm)
  
  
  message("  Ecospat")
  
  message("  ... PCA")
  
  ## PCA scores for parent1
  scores.parent1 <- suprow(pca.env, parent1_env)$li
  ## PCA scores for hybrid
  scores.hybrid <- suprow(pca.env, hybrid_env)$li
  
  message("  ... env.range.area")
  #Environmental ranges
  env_range_area <- env.range.area(scores.hybrid)
  env_range_area_P1 <- env.range.area(scores.parent1)
  
  message("  ... grid.clim.dyn")
  #Environmental overlap
  grid.clim.parent1 <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                          glob1 = scores.globclim,
                                          sp = scores.parent1, R=100,
                                          th.sp=0)
  
  
  grid.clim.hybrid <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                           glob1 = scores.globclim,
                                           sp = scores.hybrid, R=100,
                                           th.sp=0)
  
  message("  ... ecospat.niche.overlap")  
  D.overlap <- ecospat.niche.overlap(grid.clim.parent1, grid.clim.hybrid,
                                     cor = TRUE)
  Schoeners_metric <- D.overlap$D
  Hellinger_metric <- D.overlap$I
  
  jpeg(filename = paste("Plots/niche_overlap_for_", hybrid_name, "+_", parent1_name,".jpeg", sep= ""))
  
  message("  ... ecospat.plot.niche.dyn")  
  ecospat.plot.niche.dyn(grid.clim.parent1, grid.clim.hybrid, quant = 0.25,
                         interest = 2, title = "Niche Overlap",
                         name.axis1 = "PC1", name.axis2 = "PC2")
  dev.off() #prevents the display of the plots 
  
  message("  ... ecospat.niche.dyn.index")  
  #determining niche overlap b/w parent 1 and hybrid
  index_parent1comp <- ecospat.niche.dyn.index(grid.clim.parent1, grid.clim.hybrid) 
  comp1_overlap <- index_parent1comp$dynamic.index.w["stability"]
  comp1_hybrid_only <- index_parent1comp$dynamic.index.w["expansion"]
  parent1_only <- index_parent1comp$dynamic.index.w["unfilling"]
  
  message("  Geographic Ecospat")  
  ## Ecospat geographic grids
  geoGrid <- expand.grid(longitude =
                           seq(-180, 180, length.out = 200),
                         latitude =
                           seq(-90, 90, length.out = 200))
  
  grid.geo.hybrid <- ecospat.grid.clim.dyn(glob = geoGrid,
                                           glob1 = geoGrid,
                                           sp = hybridC, R=200,
                                           th.sp=0, geomask = wrld_simpl)
  
  grid.geo.parent1 <- ecospat.grid.clim.dyn(glob = geoGrid,
                                            glob1 = geoGrid,
                                            sp = parent1C, R=200,
                                            th.sp=0, geomask = wrld_simpl) 
  
  D.geo.overlap <- ecospat.niche.overlap(grid.geo.parent1, grid.geo.hybrid,
                                         cor = TRUE)
  Schoeners_metricGeo <- D.geo.overlap$D
  Hellinger_metricGeo <- D.geo.overlap$I
  
  jpeg(filename = paste("Plots/geog_", hybrid_name, "+_",
                        parent1_name,".jpeg", sep= "")) 
  
  ecospat.plot.niche.dyn(grid.geo.parent1, grid.geo.hybrid, quant = 0,
                         interest = 2, title = "Niche Overlap",
                         name.axis1 = "PC1", name.axis2 = "PC2")
  dev.off() #prevents the display of the plots 
  
  ## determining niche overlap b/w parent 1 and hybrid
  index_parent1geo <- ecospat.niche.dyn.index(grid.geo.parent1, grid.geo.hybrid)
  comp1_overlapGeo <- index_parent1geo$dynamic.index.w["stability"]
  comp1_hybrid_onlyGeo <- index_parent1geo$dynamic.index.w["expansion"]
  parent1_onlyGeo <- index_parent1geo$dynamic.index.w["unfilling"]
  
  #If the hybrid has a parent 2 then redo the above calculations with parent 2:
  if (!is.na(parent2s[i])) {
    message("  Parent 2")  
    parent2_name <- parent2s[i]
    parent2 <- read.csv(paste("Datasets/", parent2s[i], ".csv", sep=""))
    parent2C <- getcoords(parent2)
    parent2C <- SpatialPoints(parent2C, proj4string=CRS(proj4string(wrld_simpl)))
    parent2C <- parent2C[wrld_simpl] #remove points not on land
    parent2C <- data.frame(coordinates(parent2C))
    
    #Determine max and min latitude 
    maxlat_parent2 <- max(parent2C$decimalLatitude)
    
    minlat_parent2 <- min(parent2C$decimalLatitude)
    
    #determine centroid for parent 2 and position for species with two parents
    parent2_centroid <- centroid(parent2C)
    P2_centroid_latitude <- parent2_centroid[2]
    P2_centroid_longitude <- parent2_centroid[1]
    position <- "NA" #remove value from last loop
    
    if (abs(hybrid_centroid[2])>abs(parent1_centroid[2])) {
      if (abs(hybrid_centroid[2])>abs(parent2_centroid[2])){
        position <- "more polar" 
      }
    } 
    
    if (abs(hybrid_centroid[2])<abs(parent1_centroid[2])){
      if (abs(hybrid_centroid[2])<abs(parent2_centroid[2])){
        position <- "less polar"
      }
    }
    
    if (position != "more polar" & position != "less polar") {
      position <- "intermediate"
    } 
  }
    
    message("..extracting env")
    ## Extract environmental occurrence conditions
    parent2_env <- envextract(parent2C)
    
    #max, min and breadths for temperature and precipitation
    parent2_temptest <- tempbreadth(parent2_env)
    parent2_tempmax <- parent2_temptest[1]
    parent2_tempmin <- parent2_temptest[2]
    parent2_tempbreadth <- parent2_temptest[3]
    
    parent2_prectest <- precipbreadth(parent2_env)
    parent2_precmax <- parent2_prectest[1]
    parent2_precmin <- parent2_prectest[2]
    parent2_precbreadth <- parent2_prectest[3]
  
    #message("..range size & overlap")    
    ## geographical range size and overlap
    ncells <- NULL
    parent2_enm <- enmtools.species(presence.points = parent2C)
    parent2_enm$range <-
      background.raster.buffer(parent2_enm$presence.points, 100000, wc) 
    
    ncells <- freq(parent2_enm$range, value = 1, useNA= "no")
    area_parent2 <- res_km[1] *res_km[2] * ncells 
    
    comp2_overlap_geo <- geog.range.overlap(hybrid_enm, parent2_enm)
    
    ## PCA scores for parent2
    scores.parent2 <- suprow(pca.env, parent2_env)$li
    
    ## Environmental range
    env_range_area_P2 <- env.range.area(scores.parent2)
    
    #message("..ecospat")        
    ## Environmental range overlap 
    grid.clim.parent2 <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                               glob1 = scores.globclim,
                                               sp = scores.parent2, R=100,
                                               th.sp=0)
    
    D.overlap.2 <- ecospat.niche.overlap(grid.clim.parent2, grid.clim.hybrid,
                                         cor = TRUE)
    Schoeners_metric_2 <- D.overlap.2$D
    Hellinger_metric_2 <- D.overlap.2$I
    
    jpeg(filename = paste("Plots/niche_", hybrid_name, "+_",
                          parent2_name,".jpeg", sep= "")) 
    
    ecospat.plot.niche.dyn(grid.clim.parent2, grid.clim.hybrid, quant = 0.25,
                           interest = 2, title = "Niche Overlap",
                           name.axis1 = "PC1", name.axis2 = "PC2")
    dev.off()
    
    
    ## determining niche overlap b/w parent 2 and hybrid
    index_parent2comp <- ecospat.niche.dyn.index(grid.clim.parent2, grid.clim.hybrid)
    comp2_overlap <- index_parent2comp$dynamic.index.w["stability"]
    comp2_hybrid_only <- index_parent2comp$dynamic.index.w["expansion"]
    parent2_only <- index_parent2comp$dynamic.index.w["unfilling"]
    
    ###between parents
    index_parentscomp <- ecospat.niche.dyn.index(grid.clim.parent1, grid.clim.parent2)
    comp3_overlap <- index_parentscomp$dynamic.index.w["stability"]
    comp3_parent2_only <- index_parentscomp$dynamic.index.w["expansion"]
    comp3_parent1_only <- index_parentscomp$dynamic.index.w["unfilling"]
    
    #message(".. geographic ecospat")        
    ## Geographic grids
    grid.geo.parent2 <- ecospat.grid.clim.dyn(glob = geoGrid,
                                              glob1 = geoGrid,
                                              sp = parent2C, R=200,
                                              th.sp=0, geomask = wrld_simpl) 
    
    D.geo.overlap.2 <- ecospat.niche.overlap(grid.geo.parent2, grid.geo.hybrid,
                                             cor = TRUE)
    Schoeners_metricGeo_2 <- D.geo.overlap.2$D
    Hellinger_metricGeo_2 <- D.geo.overlap.2$I
    
    jpeg(filename = paste("Plots/geog_", hybrid_name, "+_",
                          parent2_name,".jpeg", sep= "")) 
    
    ecospat.plot.niche.dyn(grid.geo.parent2, grid.geo.hybrid, quant = 0,
                           interest = 2, title = "Niche Overlap",
                           name.axis1 = "PC1", name.axis2 = "PC2")
    dev.off() #prevents the display of the plots 
    
    ## determining niche overlap b/w parent 1 and hybrid
    index_parent2geo <- ecospat.niche.dyn.index(grid.geo.parent2, grid.geo.hybrid) 
    comp2_overlapGeo <- index_parent2geo$dynamic.index.w["stability"]
    comp2_hybrid_onlyGeo <- index_parent2geo$dynamic.index.w["expansion"]
    parent2_onlyGeo <- index_parent2geo$dynamic.index.w["unfilling"]
    
    ###between parents 
    index_parentsgeo <- ecospat.niche.dyn.index(grid.geo.parent1, grid.geo.parent2) 
    comp3_overlapGeo <- index_parentsgeo$dynamic.index.w["stability"]
    comp3_parent2_onlyGeo <- index_parentsgeo$dynamic.index.w["expansion"]
    comp3_parent1_onlyGeo <- index_parentsgeo$dynamic.index.w["unfilling"]
  }
  
  
  #If the hybrid does not have a parent 2 then define all of the variables as "NA"
  if (is.na(parent2s[i])){
    parent2_name <- NA
    minlat_parent2 <- NA
    maxlat_parent2 <- NA
    P2_centroid_longitude <- NA
    P2_centroid_latitude <- NA
    area_parent2 <- NA
    env_range_area_P2 <- NA
    comp2_overlap <- NA
    parent1_only <- NA
    comp2_hybrid_only <- NA
    Schoeners_metric_2 <- NA
    Hellinger_metric_2 <- NA
    comp2_overlap_geo <- NA
    parent2_temptest <- NA
    parent2_tempmax <- NA
    parent2_tempmin <- NA
    parent2_tempbreadth <- NA
    parent2_prectest <- NA
    parent2_precmax <- NA
    parent2_precmin <- NA
    parent2_precbreadth <- NA
    Schoeners_metricGeo_2 <- NA
    Hellinger_metricGeo_2 <- NA
    comp3_overlap <- NA
    parent2_only <- NA
    comp3_parent1_only <- NA
    comp2_overlapGeo <- NA
    comp2_hybrid_onlyGeo <- NA
    parent2_onlyGeo <- NA
    comp3_overlap <- NA
    comp3_parent2_only <- NA
    comp3_parent1_only <- NA
    comp3_overlapGeo <- NA
    comp3_parent2_onlyGeo <- NA
    comp3_parent1_onlyGeo <- NA
    
    if (is.na(parent2s[i])){
      parent2_name <- NA
      P2_centroid_longitude <- NA
      P2_centroid_latitude <- NA
      position <- "NA"
      
      #determine position for species with one parent
      if (abs(hybrid_centroid[2])>abs(parent1_centroid[2])) {
        position <- "more polar" 
      } else {
        position <- "less polar"
      }
      
  }
  
  data.line <- data.frame(Species = hybrid_name, P1 = parent1_name, P2 = parent2_name, 
                          Max_Lat_hybrid = maxlat_hybrid, Min_Lat_hybrid = minlat_hybrid,
                          Max_Lat_P1 = maxlat_parent1, Min_Lat_P1 = minlat_parent1,
                          Max_Lat_P2 = maxlat_parent2, Min_Lat_P2 = minlat_parent2,
                          Temp_max = hybrid_tempmax, P1_Temp_max = parent1_tempmax,
                          P2_Temp_max = parent2_tempmax, 
                          Temp_min= hybrid_tempmin, P1_temp_min = parent1_tempmin,
                          P2_temp_min = parent2_tempmin,
                          Temp_breadth = hybrid_tempbreadth, P1_temp_breadth = parent1_tempbreadth,
                          P2_temp_breadth = parent2_tempbreadth,
                          Precip_max = hybrid_precmax, P1_precip_max = parent1_precmax,
                          P2_precip_max = parent2_precmax,
                          Precip_min = hybrid_precmin, P1_precip_min = parent1_precmin,
                          P2_precip_min = parent2_precmin,
                          Precip_breadth = hybrid_precbreadth, P1_Precip_breadth = parent1_precbreadth,
                          P2_precip_breadth = parent2_precbreadth,
                          Range_km2 = area_hybrid, P1_range_km2 = area_parent1,
                          P1_hybrid_overlap_geo = comp1_overlap_geo,
                          P2_range_km2 = area_parent2, P2_hybrid_overlap_geo = comp2_overlap_geo,
                          Env_range_area = env_range_area,
                          Env_range_area_P1 = env_range_area_P1,
                          Schoeners_overlap_metric_P1H = Schoeners_metric,
                          Hellinger_metric_P1H = Hellinger_metric,
                          Env_range_area_P2 = env_range_area_P2,
                          Schoeners_overlap_metric_P2H = Schoeners_metric_2,
                          Hellinger_metric_P2H = Hellinger_metric_2,
                          P1_hybrid_overlap = comp1_overlap, 
                          Hybrid_only1 = comp1_hybrid_only, P1_only = parent1_only, 
                          P2_hybrid_overlap = comp2_overlap,
                          Hybrid_only2 = comp2_hybrid_only, P2_only = parent2_only,
                          Schoeners_overlap_geo_P1H = Schoeners_metricGeo,
                          Hellinger_geo_P1H = Hellinger_metricGeo,
                          P1_h_geo_overlap = comp1_overlapGeo,
                          Hybrid_only_geo1 = comp1_hybrid_onlyGeo,
                          P1_only_geo = parent1_onlyGeo,
                          Schoeners_overlap_geo_P2H = Schoeners_metricGeo_2,
                          Hellinger_geo_P2H = Hellinger_metricGeo_2,
                          P2_h_geo_overlap = comp2_overlapGeo,
                          Hybrid_only_geo2 = comp2_hybrid_onlyGeo,
                          P2_only_geo = parent2_onlyGeo,
                          Parents_env_overlap = comp3_overlap,
                          Parents_env_parent2only = comp3_parent2_only,
                          Parents_env_parent1only = comp3_parent1_only,
                          Parents_geo_overlap = comp3_overlapGeo,
                          Parents_geo_overlap_parent2_only = comp3_parent2_onlyGeo,
                          Parents_geo_overlap_parent1_only = comp3_parent1_onlyGeo, 
                          centroid_longitude, centroid_latitude, P1_centroid_longitude, 
                          P1_centroid_latitude, P2_centroid_longitude, P2_centroid_latitude,
                          position)
  
  
  Summary.File <- rbind(Summary.File, data.line)
}

row.names(Summary.File) <- Summary.File$Species
write.csv(Summary.File, row.names = F, file = 
            (paste("C:/Users/julia/Documents/Agriculture and Agri-Food Canada/R/polyploidy_project/", 
                   "Ploidy_project_analysis_", date.str, ".csv", sep="")))



