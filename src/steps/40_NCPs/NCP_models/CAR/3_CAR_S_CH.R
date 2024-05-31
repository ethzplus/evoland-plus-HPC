#####
#
# The aim of this script is to reunite the different InVEST models for each
# region previously calculated in the other code files 1_CAR... & 2_CAR_...
# into one final map
#
#####

library(terra)

# Load the parameters into env by sourcing the ../load_params.R script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.dir <- dirname(sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]))
source(file.path(script.dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$CAR$out_prefix)) { stop("params$CAR$out_prefix is not set") }
if (is.null(params$run_params$NCP_RUN_SCENARIO_ID)) {
  stop("params$run_params$NCP_RUN_SCENARIO_ID is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$NCP_RUN_OUTPUT_DIR)) {
  stop("params$run_params$NCP_RUN_OUTPUT_DIR is not set")
}
if (is.null(params$run_params$NCP_RUN_SCRATCH_DIR)) {
  stop("params$run_params$NCP_RUN_SCRATCH_DIR is not set")
}
if (is.null(params$data$lulc)) { stop("params$data$lulc is not set") }


# inv_model_dir <- params$CAR$dir18 # folder containing invest models
inv_model_dir <- file.path(params$run_params$NCP_RUN_SCRATCH_DIR,
                           params$run_params$NCP_RUN_SCENARIO_ID, "CAR",
                           params$run_params$NCP_RUN_YEAR, "Invest_models")
# results <- params$CAR$res18 # folder to save results
results <- file.path(params$run_params$NCP_RUN_OUTPUT_DIR,
                     params$run_params$NCP_RUN_SCENARIO_ID, "CAR",
                     params$run_params$NCP_RUN_YEAR)
# n18 <- params$CAR$n18 # name of the final raster
name_out <- paste(params$CAR$out_prefix, params$run_params$NCP_RUN_YEAR,
                  sep = "")

# Function to apply for each time period

carbon_process_3 <- function(ind_dir, raster_out) {

  # for each folder in list_dirs, search file matching 'tot_c_cur_*.tif'
  list_files <- list.files(ind_dir, full.names = T, recursive = T, pattern = "tot_c_cur_")
  # exclude all files inside tot_c_united
  list_files <- list_files[!grepl("tot_c_united", list_files)]

  #--- bind together each file
  # creating empty raster to fill
  # TODO: Change for different countries
  bind <- rast()
  crs(bind) <- params$proj$crs
  ext(bind) <- c(params$proj$ext)
  res(bind) <- params$proj$res

  print(paste("Using resolution ", res(bind), " for final raster"))

  # merging each individual carbon stock raster into the empty raster created
  # above

  for (i in 1:length(list_files)) {
    name <- list_files[i]
    # loading a raster file from the folder
    region_rast <- rast(name)
    crs(region_rast) <- params$proj$crs
    # check if resolution is equal to the one of the empty raster
    if (any(res(region_rast) != res(bind))) {
      print(paste("Resolution of ", name, " is not equal to the final raster"))
      print(paste("Resolution of ", name, " is ", res(region_rast)))
      print(paste("Resolution of final raster is ", res(bind)))
      print("Exiting")
      stop("Resolution of raster is not equal to the final raster")
    }
    # binding using the max function -> if there is an overlap, the max value is preserved
    bind <- mosaic(bind, region_rast, fun = "max")
    print(paste(name, " added to final raster"))
  }

  crs(bind) <- params$proj$crs
  bind <- extend(bind, c(params$proj$ext))  # extend to the extent of the final raster

  writeRaster(bind, raster_out, overwrite = TRUE)

  cat(paste("Final raster saved to ", raster_out, "\n"))
}

#--- Applying for each period

carbon_process_3(
  inv_model_dir,
  file.path(results, paste(name_out, ".tif", sep = ""))
)
