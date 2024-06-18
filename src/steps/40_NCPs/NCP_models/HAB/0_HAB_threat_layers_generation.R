#########################################################################
## --- Creating threat raster layers for Habitat quality module InVEST
## --- These are based on LULC maps and municipalities polygons

## The Swiss municipality layer is processed to calculate population density
## and identify rural residential areas.
## The threat_hab function is then defined to create the threat layers.
## This function classifies the land use/land cover raster into crop, rural
## residential, and urban threat layers.
#########################################################################

library(raster)
library(terra)
library(rgdal)
library(sf)

# Load the parameters into env by sourcing the ../load_params.R script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_dir <- dirname(sub(
  file_arg_name, "",
  initial_options[grep(file_arg_name, initial_options)]
))
source(file.path(script_dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$run_params$NCP_RUN_SCENARIO_ID)) {
  stop("params$run_params$NCP_RUN_SCENARIO_ID is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$NCP_RUN_SCRATCH_DIR)) {
  stop("params$run_params$NCP_RUN_SCRATCH_DIR is not set")
}
if (is.null(params$data$lulc)) {
  stop("params$data$lulc is not set")
}
if (is.null(params$data$rur_res)) {
  stop("params$data$rur_res is not set")
}
if (is.null(params$data$urb_res)) {
  stop("params$data$urb_res is not set")
}

#--Local variables
out_fold <- file.path(
  params$run_params$NCP_RUN_SCRATCH_DIR,
  params$run_params$NCP_RUN_SCENARIO_ID, "HAB",
  params$run_params$NCP_RUN_YEAR
)

# LULC raster
lulc <- rast(params$data$lulc)
rur_res <- rast(params$data$rur_res)
urb_res <- rast(params$data$urb_res)

# Function to generate threat layers (crops, urban and rural residential)
threat_hab <- function(lulc, rur_res, urb_res, out_fold) {
  timestart <- Sys.time()
  #-Crop threat layer
  # reclassification scheme using a range e.g. 36-41 and new value 1
  # TODO: needs change to use the lulcc classification
  #       upper value in range include, lower values is not
  c1 <- c(
    0, 14, 0,
    14, 15, 1,
    15, 17, 0,
    17, 18, 1,
    18, 21, 0
  )
  mc1 <- matrix(c1, ncol = 3, byrow = TRUE)
  crop <- terra::classify(lulc, mc1) # nolint

  # Export
  terra::writeRaster(crop, file.path(out_fold, "crop_c.tif"), overwrite = TRUE)
  print("Crop layer created")

  ### -Rural residential

  # mask lulc so that only lulc values within rural areas are retained
  rur_lulc <- terra::mask(lulc, rur_res)

  #-selecting residential areas
  # reclassification scheme for residential areas within rural areas
  r1 <- c(
    0, 9, 0,
    9, 10, 1,
    10, 21, 0
  )

  mr1 <- matrix(r1, ncol = 3, byrow = TRUE)
  rures <- terra::classify(rur_lulc, mr1)
  # filling gaps with 0 on entire extent
  exp_rr <- terra::classify(rures, cbind(NA, 0))

  # Export
  terra::writeRaster(
    exp_rr, file.path(out_fold, "rures_c.tif"),
    overwrite = TRUE
  )
  print("Rural residential layer created")

  ### -urban

  # mask so that only values of raster cells within urban residential are
  # retained
  urb_lulc <- terra::mask(lulc, urb_res)

  # reclassify using same scheme as for rural residential
  urbres <- terra::classify(urb_lulc, mr1)
  # filling gaps with 0 on all extent
  exp_ur <- terra::classify(urbres, cbind(NA, 0))

  # Export
  terra::writeRaster(
    exp_ur, file.path(out_fold, "urban_c.tif"),
    overwrite = TRUE
  )
  print("Urban layer created")
  timestop <- Sys.time()
  elapsed <- difftime(timestop, timestart, units = "secs")
  print(paste("duration: ", elapsed[[1]], " seconds"))
}

# Create the output folder if it does not exist
if (!dir.exists(out_fold)) {
  dir.create(out_fold, recursive = TRUE)
}
#### Run the function on the input data
threat_hab(lulc, rur_res, urb_res, out_fold)
