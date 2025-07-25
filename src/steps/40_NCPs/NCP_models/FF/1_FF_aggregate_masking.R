###################
#
# The aim of this script is to
# 1) aggregate map outputs from ecocrop models,
# 2) mask aggregated map with only agricultural land
#
# The code processes a land use/land cover raster to identify agricultural
# areas and then uses this information to mask a stack of Ecocrop map outputs.
# The masked raster values are then normalized to a range between 0 and 1.
###################

library(terra)

# Load the parameters into env by sourcing the ../load_params.R script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_dir <- dirname(sub(
  file_arg_name, "", initial_options[grep(file_arg_name, initial_options)]
))
source(file.path(script_dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$proj$crs)) {
  stop("params$proj$crs is not set")
}
if (is.null(params$FF$ecocrop_dir)) {
  stop("params$FF$ecocrop_dir is not set")
}
if (is.null(params$run_params$NCP_RUN_SCENARIO_ID)) {
  stop("params$run_params$NCP_RUN_SCENARIO_ID is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$NCP_RUN_RCP)) {
  stop("params$run_params$NCP_RUN_RCP is not set")
}
if (is.null(params$run_params$NCP_RUN_OUTPUT_DIR)) {
  stop("params$run_params$NCP_RUN_OUTPUT_DIR is not set")
}
if (is.null(params$data$lulc)) {
  stop("params$data$lulc is not set")
}

# Local variables
# Land use/land cover map
lulc <- rast(params$data$lulc)
# Ecocrop map outputs
eco_maps <- rast(list.files(
  file.path(
    params$FF$ecocrop_dir, params$run_params$NCP_RUN_YEAR,
    params$run_params$NCP_RUN_RCP
  ),
  full.names = TRUE, pattern = "\\.tif$"
))
results <- file.path(
  params$run_params$NCP_RUN_OUTPUT_DIR,
  params$run_params$NCP_RUN_SCENARIO_ID, "FF"
)
dir.create(results, showWarnings = FALSE, recursive = TRUE)
name_out <- paste("FF_S_CH_", params$run_params$NCP_RUN_YEAR, ".tif", sep = "")

# Reclassify agriculture categories of landuse to 1, rest to 0
# the first two numbers describe a range and the third the new value
m <- c(0, 15, 0,
       15, 19, 1,
       19, Inf, 0)

m1 <- matrix(m, byrow = TRUE, ncol = 3)
lulc_agr <- classify(lulc, m1)

# adapt extent
lulc_agr <- crop(lulc_agr, ext(eco_maps[[1]]))
lulc_agr <- extend(lulc_agr, ext(eco_maps[[1]]))
crs(lulc_agr) <- params$proj$crs

# Mask for agricultural areas and apply a mean
ff_agr <- mean(eco_maps)
ff_agr_masked <- lulc_agr * ff_agr


# Normalize value
nx <- minmax(ff_agr_masked)
rn <- (ff_agr_masked - nx[1, ]) / (nx[2, ] - nx[1, ])

# Export
writeRaster(rn, file.path(results, name_out), overwrite = TRUE)
print(paste("FF supply layer created: ", file.path(results, name_out)))
