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

# Load the parameters into env by sourcing the ../load_params.R script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_dir <- dirname(sub(file_arg_name, "", initial_options[
  grep(file_arg_name, initial_options)
]))
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
if (is.null(params$run_params$NCP_RUN_OUTPUT_DIR)) {
  stop("params$run_params$NCP_RUN_OUTPUT_DIR is not set")
}

data_folder <- file.path(
  params$run_params$NCP_RUN_OUTPUT_DIR,
  params$run_params$NCP_RUN_SCENARIO_ID,
  "POL",
  params$run_params$NCP_RUN_YEAR
)

# Retrieve a list of file paths from the data folder that match the pattern
files <- list.files(data_folder, pattern = "supply", full.names = TRUE)

out_dir <- file.path(
  params$run_params$NCP_RUN_OUTPUT_DIR,
  params$run_params$NCP_RUN_SCENARIO_ID,
  "POL"
)
# Create the output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
name_out <- paste("POL_S_CH_", params$run_params$NCP_RUN_YEAR, ".tif", sep = "")

# Read and stack the rasters from the list of file paths into a single raststack
p_list <- rast(files)

# Calculate the sum of the values across all rasters in the stack
pol_sum <- sum(p_list)

# Normalize the summed raster values to a range of 0 to 1
nx <- minmax(pol_sum)
rn <- (pol_sum - nx[1,]) / (nx[2,] - nx[1,]) # nolint

terra::writeRaster(rn, file.path(out_dir, name_out), overwrite = TRUE)
print(paste("Pollination supply layer created: ", file.path(out_dir, name_out)))

# Remove temporary files not explicitly needed,
# as they are stored in the scratch directory
