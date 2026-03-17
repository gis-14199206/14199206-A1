# Course Number: GEOG71922; Student ID Number:14199206

# Outline:
# 1. Read in and map the data
# 2. Create a vector layer from coordinates
# 3. Create spatial objects and define study extent
# 4. Load and prepare environmental rasters
# 5. Reclassify land cover and derive habitat variables
# 6. Generate background points
# 7. Conduct scale analysis for broadleaf woodland and select optimum neighbourhood size
# 8. Build final environmental predictor stack
# 9. Extract covariates to presence and background points
# 10. Explore predictor relationships and inspect data structure
# 11. Fit GLM
# 12. Fit Maxnet
# 13. Validate models using cross-validation
# 14. Compare model performance
# 15. Predict habitat suitability across the study area
# 16. Plot final maps and interpretation figures



# 1. Read in and map the data

# set working directory
setwd("D:/71922")

#Install packages
install.packages(c("terra","sf","mapview"))

# load libraries
library(terra)
library(sf)
library(dismo)
library(maxnet)
library(glmnet)
library(precrec)
library(ggplot2)

# read in the spreadsheet containing the coordinates for occurrence locations
meles <- read.csv("Melesmeles.csv")

# view the first six rows of the data
head(meles)

# check the class of the loaded .csv file
class(meles)

# subset the data to only include points with complete coordinates
meles <- meles[!is.na(meles$Latitude), ]

# remove all records with coordinate uncertainty greater than 1000 m
meles <- meles[meles$Coordinate.uncertainty_m <= 1000, ]

# check number of rows
nrow(meles)





# 2. Create a vector layer from coordinates

# make spatial points layer
# create coordinates object
meles.latlong <- data.frame(x = meles$Longitude, y = meles$Latitude)

# use coordinates object to create spatial points object
meles.sp <- vect(meles.latlong, geom = c("x","y"))

# check that the points now have our desired CRS
crs(meles.sp) <- "epsg:4326"

# plot our points
plot(meles.sp)

# 3. Create spatial objects and define study extent（week3）

# read in the study area polygon
scot <- st_read("scotSamp.shp")

# inspect study area geometry only
plot(st_geometry(scot))

# project badger points to the same CRS as the study area
meles.sp <- project(meles.sp, st_crs(scot)$wkt)

# crop badger points to the study area
meles.fin <- crop(meles.sp, vect(scot))


# inspect
plot(st_geometry(scot))
plot(meles.fin, add = TRUE)

# 4. Load and prepare environmental rasters

# load in the land cover map
LCM <- rast("LCMUK.tif")

# select the land cover classification layer
LCM <- LCM$LCMUK_1

# crop to the extent of the study area plus a little more
LCM <- crop(LCM, vect(st_buffer(scot, dist = 1000)))

# now mask to the actual study area boundary
LCM <- mask(LCM, vect(scot))

# inspect
plot(LCM)
plot(meles.fin, add = TRUE)

# load elevation data (Week 3 practical)
DEM <- rast("demScotland.tif")

# crop to study area
DEM <- crop(DEM, vect(st_buffer(scot, dist = 1000)))

# mask to study area boundary
DEM <- mask(DEM, vect(scot))

# inspect
plot(DEM)
plot(meles.fin, add = TRUE)

# 5. Reclassify land cover and derive habitat variables

# treat land cover raster as categorical
LCM <- as.factor(LCM)

# inspect land cover classes
levels(LCM)

# -------------------------
# 5a. Broadleaf woodland
# -------------------------

# create reclassification vector:
# broadleaf woodland = class 1
# all other classes = 0
reclass <- c(0, 1, rep(0, nrow(levels(LCM)[[1]]) - 2))

# combine with original classes
RCmatrix <- cbind(levels(LCM)[[1]], reclass)

# keep only class code + new value
RCmatrix <- RCmatrix[, 2:3]

# convert to numeric
RCmatrix <- apply(RCmatrix, 2, FUN = as.numeric)

# reclassify raster
broadleaf <- classify(LCM, RCmatrix)

# inspect
plot(broadleaf, main = "Broadleaf woodland")
plot(meles.fin, add = TRUE, col = "red")


# -------------------------
# 5b. Urban land cover
# -------------------------

# create reclassification vector:
# urban + suburban = 1 (assumed last two classes)
reclassUrban <- c(rep(0, nrow(levels(LCM)[[1]]) - 2), 1, 1)

# combine with original classes
RCmatrixUrban <- cbind(levels(LCM)[[1]], reclassUrban)

# keep only class code + new value
RCmatrixUrban <- RCmatrixUrban[, 2:3]

# convert to numeric
RCmatrixUrban <- apply(RCmatrixUrban, 2, FUN = as.numeric)

# reclassify raster
urban <- classify(LCM, RCmatrixUrban)

# inspect
plot(urban, main = "Urban land cover")
plot(meles.fin, add = TRUE, col = "blue")


# 6. Generate background points

# set seed for reproducibility
set.seed(11)

# sample random background points from the study area
back.xy <- spatSample(broadleaf, size = 1000, as.points = TRUE, na.rm = TRUE)

# inspect
plot(broadleaf)
plot(meles.fin, add = TRUE, col = "red")
plot(back.xy, add = TRUE, col = "blue")

# 7. Extract covariates to points

# extract broadleaf values to background points
eA <- extract(broadleaf, back.xy)

# extract broadleaf values to presence points
eP <- extract(broadleaf, meles.fin)

# inspect frequency tables
table(eA[,2])
table(eP[,2])

# create data frames for absence and presence
Abs <- data.frame(crds(back.xy), broadleaf = eA[,2], Pres = 0)
Pres <- data.frame(crds(meles.fin), broadleaf = eP[,2], Pres = 1)

# combine presence and background data
melesData <- rbind(Pres, Abs)

# inspect
head(melesData)
summary(melesData)

# 8. Fit GLM

# inspect the response and predictor
table(melesData$Pres)
table(melesData$broadleaf)
table(melesData$broadleaf, melesData$Pres)

# fit a binomial GLM
glm.meles <- glm(Pres ~ broadleaf, family = "binomial", data = melesData)

# inspect model summary
summary(glm.meles)

# inspect model fit using AIC
AIC(glm.meles)

# generate predicted probabilities for the point data
melesData$glm.pred <- predict(glm.meles, type = "response")

# inspect
head(melesData)
summary(melesData$glm.pred)

# 9. Fit Maxnet model

library(maxnet)

# response variable
p <- melesData$Pres

# predictor
env <- data.frame(broadleaf = melesData$broadleaf)

# fit model
maxnet.meles <- maxnet(p = p, data = env, f = maxnet.formula(p, env))

# check model
maxnet.meles

# predict probability
melesData$maxnet.pred <- predict(maxnet.meles, env, type = "cloglog")

# inspect
summary(melesData$maxnet.pred)
