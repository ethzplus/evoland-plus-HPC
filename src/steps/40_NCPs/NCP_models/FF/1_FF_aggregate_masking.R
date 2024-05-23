###################
#
# The aim of this script is to 1) aggregate map outputs from ecocrop models, 2) mask aggregated map with only agricultural land
#
#The code processes a land use/land cover raster to identify agricultural areas and then uses this information to mask a stack of Ecocrop map outputs. 
#The masked raster values are then normalized to a range between 0 and 1.
###################

library(terra)

# Load the parameters into env by sourcing the ../load_params.R script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.dir <- dirname(sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]))
source(file.path(script.dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$proj$crs)) { stop("params$proj$crs is not set") }
if (is.null(params$FF$ecocrop_dir)) { stop("params$FF$ecocrop_dir is not set") }
if (is.null(params$run_params$NCP_RUN_LULC_MAP)) {
  stop("params$run_params$NCP_RUN_LULC_MAP is not set")
}
if (is.null(params$run_params$NCP_RUN_OUTPUT_DIR)) {
  stop("params$run_params$NCP_RUN_OUTPUT_DIR is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}

# Local variables
# Land use/land cover map
lulc <- rast(params$run_params$NCP_RUN_LULC_MAP)
# Ecocrop map outputs
eco_maps <- rast(list.files(
  file.path(params$FF$ecocrop_dir), full.names = T, pattern = "\\.tif$"
))
results <- file.path(params$run_params$NCP_RUN_OUTPUT_DIR, "FF",
                     params$run_params$NCP_RUN_YEAR)
dir.create(results, showWarnings = FALSE, recursive = TRUE)

# Reclassify agriculture categories of landuse to 1, rest to 0
# TODO: needs change to use the LULCC classes
# in the example below the first two numbers describe a range and the third the new value
m <- c(0, 37, 0,
       37, 48, 1,
       48, Inf, 0)

m1 <- matrix(m, byrow = TRUE, ncol = 3)
lulc_agr <- classify(lulc, m1)

# adapt extent
lulc_agr <- crop(lulc_agr, ext(eco_maps[[1]]))
lulc_agr <- extend(lulc_agr, ext(eco_maps[[1]]))
crs(lulc_agr) <- crs(params$proj$crs)

# Mask for agricultural areas and apply a mean
ff_agr <- mean(eco_maps)
ff_agr_masked <- lulc_agr * ff_agr


# Normalize value
nx <- minmax(ff_agr_masked)
rn <- (ff_agr_masked - nx[1,]) / (nx[2,] - nx[1,])

#Export
writeRaster(rn, file.path(results, "FF_S_CH.tif"), overwrite = TRUE)
print(paste("Exported FF_S_CH.tif to", results))
