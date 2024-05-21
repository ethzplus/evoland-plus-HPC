################
#
# Code to put together pollination information
#
# The code stacks all files of the data_folder with pattern = supply and then
# calculates the sum and mean of the raster stack, aggregating the pollination
# information. Finally a normalization to values between 0 and 1 is conducted.
#
################


library(terra)
library(raster)

# Load the parameters into env by sourcing the ../load_params.R script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.dir <- dirname(sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]))
source(file.path(script.dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$NCP_RUN_SCRATCH_DIR)) {
  stop("params$run_params$NCP_RUN_SCRATCH_DIR is not set")
}
if (is.null(params$run_params$NCP_RUN_OUTPUT_DIR)) {
  stop("params$run_params$NCP_RUN_OUTPUT_DIR is not set")
}

data_folder <- file.path(params$run_params$NCP_RUN_SCRATCH_DIR, "POL",
                         params$run_params$NCP_RUN_YEAR)  # TODO: Check if files are here

# Retrieve a list of file paths from the data folder that match the pattern "supply"
files <- list.files(data_folder, pattern = "supply", full.names = T)

out_dir <- file.path(params$run_params$NCP_RUN_OUTPUT_DIR, "POL",
                     params$run_params$NCP_RUN_YEAR)

# Read and stack the rasters from the list of file paths into a single raster stack
p_list <- rast(files)

# Calculate the sum of the values across all rasters in the stack
pol_sum <- sum(p_list)

# Normalize the summed raster values to a range of 0 to 1
nx <- minmax(pol_sum)
rn <- (pol_sum - nx[1,]) / (nx[2,] - nx[1,])

writeRaster(rn, file.path(out_dir, "POL_S_CH.tif"))

# Remove temporary files not explicitly needed,
# as they are stored in the scratch directory