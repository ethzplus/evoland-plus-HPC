###################
#
# The aim of this script is to 1) aggregate map outputs from ecocrop models, 2) mask aggregated map with only agricultural land
#
#The code processes a land use/land cover raster to identify agricultural areas and then uses this information to mask a stack of Ecocrop map outputs. 
#The masked raster values are then normalized to a range between 0 and 1.
###################

library(terra)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

# Local variables
# Land use/land cover map
lulc18 <- rast(params$data$lulcc2018)
# Ecocrop map outputs
r_stack_2_dir <- rast(
  list.files(params$FF$r_stack_2_dir, full.names = T, pattern = "\\.tif$")
)
results <- paste(wd, "results", sep = "/") # result folder
crs(r_stack_2_dir) <- params$proj$crs

# Reclassify agriculture categories of landuse to 1, rest to 0
# TODO: needs change to use the LULCC classes
# in the example below the first two numbers describe a range and the third the new value
m <- c(0, 37, 0,
       37, 48, 1,
       48, Inf, 0)

m1 <- matrix(m, byrow = TRUE, ncol = 3)

lulc_agr <- classify(lulc18, m1)

# adapt extent

lulc_agr <- crop(lulc_agr, ext(r_stack_2_dir[[1]]))
lulc_agr <- extend(lulc_agr, ext(r_stack_2_dir[[1]]))

# Mask for agricultural areas and apply a mean

ff_agr <- mean(r_stack_2_dir)

ff_agr_masked <- lulc_agr * ff_agr


# Normalize value

nx <- minmax(ff_agr_masked)
rn <- (ff_agr_masked - nx[1,]) / (nx[2,] - nx[1,])

rn

#Export

writeRaster(rn, paste(results, "FF_S_CH.tif", sep = "/"))


