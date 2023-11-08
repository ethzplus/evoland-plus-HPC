#########################################################################
##--- Creating threat raster layers for Habitat quality module InVEST
##--- These are based on LULC maps and municipalities polygons

##. The Swiss municipality layer is processed to calculate population density and identify rural residential areas. 
##The threat_hab function is then defined to create the threat layers. 
##This function classifies the land use/land cover raster into crop, rural residential, and urban threat layers.
#########################################################################

library(raster)
library(terra)
library(rgdal)
library(sf)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

#--Local variables
setwd(params$HAB$wd)
out_fold <- paste(params$HAB$wd, "LAYERPROD", sep = "/")

# LULC raster
lulc18 <- rast(params$data$lulcc2018)
# municipality layer
mun_layer <- vect(params$data$mun_layer)


#### 1) processing of the municipality layer -- adding a population density column to the municpalities
# TODO: dependent on the municipality layer, needs to be adapted for other countries
mun_layer$EINWOHNERZ <- as.numeric(mun_layer$EINWOHNERZ)              #nb inhabitants
mun_layer$density <- mun_layer$EINWOHNERZ / (mun_layer$GEM_FLAECH / 100)  #surface area, giving density
mun_layer <- mun_layer[-which(is.na(mun_layer$density)),]             #removing NA values

#-- making the distinction between rural residential and urban residential (rural: <10000 habitants OR density <100)
# TODO: dependent on the municipality layer, needs to be adapted for other countries
rur_res <- mun_layer[mun_layer$EINWOHNERZ < 10000 | mun_layer$density < 100,] # For rural residential

#### 2) Function to generate threat layers (crops, urban and rural residential)

threat_hab <- function(lulc, year) {

  timestart <- Sys.time()
  #-Crop threat layer
  # reclassification scheme using a range e.g. 36-41 and new value 1
  # TODO: needs change to use the lulcc classification
  #       upper value in range include, lower values is not
  c1 <- c(0, 14, 0,
          14, 15, 1,
          15, 17, 0,
          17, 18, 1,
          18, 21, 0)
  mc1 <- matrix(c1, ncol = 3, byrow = TRUE)
  crop <- classify(lulc, mc1)

  cropname <- paste(year, "_crop_c.tif", sep = "")
  writeRaster(crop, paste(out_fold, cropname, sep = "/"), overwrite = TRUE) # Export
  print(paste(cropname, " layer created", sep = ""))

  ###-Rural residential

  #dissolving polygon of rural areas
  rur <- aggregate(rur_res)
  #cropping lulc to rural areas
  rur_lu <- crop(lulc, rur)
  #mask so that only lulc values within rural areas are retained
  rur_lulc <- mask(rur_lu, rur)

  #-selecting residential areas
  #reclassification scheme for residential areas within rural areas
  # TODO: needs change to use the lulcc classification
  r1 <- c(0, 9, 0,
          9, 10, 1,
          10, 21, 0)

  mr1 <- matrix(r1, ncol = 3, byrow = TRUE)
  rures <- classify(rur_lulc, mr1)
  exp_rr <- classify(rures, cbind(NA, 0)) #filling gaps with 0 on entire extent

  rrname <- paste(year, "_rures_c.tif", sep = "")
  writeRaster(exp_rr, paste(out_fold, rrname, sep = "/"), overwrite = TRUE) # Export
  print(paste(rrname, " layer created", sep = ""))
  ###-urban

  outline <- aggregate(mun_layer) #getting Switzerland blank canvas
  urb_res <- erase(outline, rur) #remove all rural residential areas, leaving us with urban residential

  # crop to urban residential
  urb_lu <- crop(lulc, urb_res)
  #mask so that only values of raster cells within urban residential are retained
  urb_lulc <- mask(rur_lu, urb_res)

  urbres <- classify(urb_lulc, mr1) #reclassify using same scheme as for rural residential
  exp_ur <- classify(urbres, cbind(NA, 0)) # filling gaps with 0 on all extent

  urname <- paste(year, "_urban_c.tif", sep = "")
  writeRaster(exp_ur, paste(out_fold, urname, sep = "/"), overwrite = TRUE) #Export
  print(paste(urname, " layer created"))
  timestop <- Sys.time()
  elapsed <- difftime(timestop, timestart, units = "mins")
  print(paste("duration: ", elapsed[[1]], " mins", sep = ""))

}


#### 3) Applying the function to desired time period. 

threat_hab(lulc18, "18")
