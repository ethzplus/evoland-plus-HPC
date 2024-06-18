###################
#
# The aim of this script is to apply an "ecocrop" model to each crop of
# interest to generate suitability maps
#
# The script reads precipitation, temperature, pH and crop datasets and
# preprocesses them.
# For each crop in the list, the script runs the Ecocrop model,
# which predicts the suitability of the crop based on the provided
# environmental variables.
###################

library(Recocrop)
library(data.table)
library(terra)
library(sf)

# Load the parameters into env by sourcing the ../load_params.R script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_dir <- dirname(sub(
  file_arg_name, "",
  initial_options[grep(file_arg_name, initial_options)]
))
source(file.path(script_dir, "..", "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$FF$crops_data)) {
  stop("params$FF$crops_data is not set")
}
if (is.null(params$FF$ecocrop_dir)) {
  stop("params$FF$ecocrop_dir is not set")
}
if (is.null(params$data$ph_raster)) {
  stop("params$data$ph_raster is not set")
}
if (is.null(params$data$pavg_dir)) {
  stop("params$data$pavg_dir is not set")
}
if (is.null(params$data$tavg_dir)) {
  stop("params$data$tavg_dir is not set")
}
if (is.null(params$proj$crs)) {
  stop("params$proj$crs is not set")
}


#--local variables
# setting paths to local variables
out_dir <- file.path(params$FF$ecocrop_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pavg <- params$data$pavg_dir # precipitations
tavg <- params$data$tavg_dir # temperature
ph <- params$data$ph_raster # pH
crops_data <- params$FF$crops_data # crops names

## ---Loading variables

# List of food of interest

# read table containing crops list

crops <- t(fread(crops_data, sep = ";", header = FALSE))
rownames(crops) <- c(seq_len(nrow(crops)))
colnames(crops) <- "crop"
crops <- as.data.frame(crops)

# Environmental variables (Temp, prec, ph)
ta <- list.files(tavg) # list of temperature files
ta <- ta[grepl(".tif$", ta)] # filtering only tif files
# TODO: Possibly needs change for other countries
# adjust temperature values if needed
# (originally multiplied by 10)
print("Listed temperature files")

pr <- list.files(pavg) # list of precipitation files
pr <- pr[grepl(".tif$", pr)] # filtering only tif files
# TODO: Possibly needs change for other countries
# adjust precipitation values if needed
# (originally multiplied by 100 and for 100m raster cells)
print("Listed precipitation files")

# check that the number of temperature and precipitation files are the same
if (length(ta) != length(pr)) {
  stop("The number of temperature and precipitation files are not the same")
}

ph <- rast(ph)
ph <- aggregate(ph, fact = 4, fun = mean)
print("Loaded and aggregated pH data")

# Loop Year & RCP combinations
# format of ta Tave_{year}_{RCPxx}_m{month}.tif
# find unique years and rcp
years <- unique(substr(ta, start = 6, stop = 9))
rcps <- unique(substr(ta, start = 11, stop = 15))

print(paste("Found", length(years), "years and", length(rcps), "rcps:"))
print(years)
print(rcps)

ecocrop_models <- function(crops, ta, pr, ph, use_crs, out_dir) {
  ## ---- Ecocrop model on each crop type
  # counter variable
  z <- 1
  for (i in seq_len(nrow(crops))) {
    # Extract the crop name from the current row.
    a <- crops$crop[i]
    # Create a short version of the name containing 3 characters
    c_name <- substr(a, start = 1, stop = 3)
    # if short version allready exists in crops$mod then z is added to the name
    if (c_name %in% crops$mod_name == TRUE) {
      crops$mod_name[i] <- paste(c_name, z, sep = "_")
      z <- z + 1
    } else {
      # mod_name becomes a column with each model name
      crops$mod_name[i] <- c_name
    }
    # Get parameters for the current crop and create a model
    pars <- Recocrop::ecocropPars(a)
    m_name <- Recocrop::ecocrop(pars) # model
    Recocrop::control(m_name, get_max = TRUE) # get the maximum value of the
    # model
    # Predict using the model
    pred <- predict(
      m_name,
      tavg = ta, prec = pr, ph = ph, wopt = list(names = c_name)
    )
    terra::crs(pred) <- terra::crs(use_crs)
    # Export the prediction
    ex_name <- paste(c_name, ".tif", sep = "")
    terra::writeRaster(pred, file.path(out_dir, ex_name), overwrite = TRUE)
    # Print progress
    print(paste(
      "model:", a, "done.", round((i * 100) / nrow(crops), 1), "%",
      sep = " "
    ))
    gc()
  }
}


# Loop through each year and RCP combination
for (year in years) {
  for (rcp in rcps) {
    # filter temperature and precipitation files for the current year and RCP
    ta_year_rcp <- ta[grepl(paste(year, rcp, sep = "_"), ta)]
    pr_year_rcp <- pr[grepl(paste(year, rcp, sep = "_"), pr)]
    print(paste(
      "Year:", year,
      "RCP:", rcp,
      "Temperature files:", length(ta_year_rcp),
      "Precipitation files:", length(pr_year_rcp)
    ))
    # check that the number of temperature and precipitation files are 12
    if (length(ta_year_rcp) != 12 || length(pr_year_rcp) != 12) {
      stop("The number of temperature and precipitation files are not 12 each")
    }
    # read raster files
    ta_year_rcp <- rast(file.path(tavg, ta_year_rcp))
    pr_year_rcp <- rast(file.path(pavg, pr_year_rcp))
    # Loop through each crop
    current_out_dir <- file.path(out_dir, year, rcp)
    dir.create(current_out_dir, showWarnings = FALSE, recursive = TRUE)
    ecocrop_models(
      crops, ta_year_rcp, pr_year_rcp, ph, params$proj$crs,
      current_out_dir
    )
  }
}
