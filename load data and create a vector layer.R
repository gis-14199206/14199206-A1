# Course Number: GEOG71922; Student ID Number:14199206

# Outline:
# 1. Read in and map the data
# 2. Create a vector layer from coordinates
# 3. Create spatial objects and define study extent
# 4. Load and prepare environmental rasters
# 5. Scale analysis for broadleaf woodland
# 6. Generate background points
# 7. Extract covariates to points
# 8. Fit GLM
# 9. Fit Maxnet / Random Forest
# 10. Cross-validation and model evaluation
# 11. Predict habitat suitability
# 12. Plot final maps



# 1. Read in and map the data

# set working directory
setwd("D:/71922")

#Install packages
install.packages(c("terra","sf","mapview"))

# load libraries
library(terra)
library(sf)

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

