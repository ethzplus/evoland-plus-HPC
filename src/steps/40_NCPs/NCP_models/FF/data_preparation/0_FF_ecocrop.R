###################
#
# The aim of this script is to apply an "ecocrop" model to each crop of interest to generate suitability maps
#
#The script reads precipitation, temperature, pH and crop datasets and preprocesses them. 
#For each crop in the list, the script runs the Ecocrop model,
#which predicts the suitability of the crop based on the provided environmental variables.
###################

library(Recocrop)
library(data.table)
library(terra)
library(sf)

# Load the parameters into env by sourcing the ../load_params.R script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.dir <- dirname(sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]))
source(file.path(script.dir, "..", "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$FF$crops_data)) { stop("params$FF$crops_data is not set") }
if (is.null(params$FF$ecocrop_dir)) { stop("params$FF$ecocrop_dir is not set") }
if (is.null(params$data$ph_raster)) { stop("params$data$ph_raster is not set") }
if (is.null(params$data$pavg_dir)) { stop("params$data$pavg_dir is not set") }
if (is.null(params$data$tavg_dir)) { stop("params$data$tavg_dir is not set") }
if (is.null(params$proj$crs)) { stop("params$proj$crs is not set") }


#--local variables
# setting paths to local variables
out_dir <- file.path(params$FF$ecocrop_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pavg95 <- params$data$pavg_dir  # precipitations
tavg95 <- params$data$tavg_dir  # temperature
ph <- params$data$ph_raster  # pH
crops_data <- params$FF$crops_data  # crops names

##---Loading variables

#List of food of interest

#read table containing crops list

crops <- t(fread(crops_data, sep = ";", header = FALSE))
rownames(crops) <- c(1:nrow(crops))
colnames(crops) <- "crop"
crops <- as.data.frame(crops)

#Environmental variables (Temp, prec, ph)
ta <- list.files(tavg95)  # list of temperature files
print("Loaded temperature data")
ta <- ta[grepl(".tif$", ta)]  # filtering only tif files
ta <- rast(paste(tavg95, ta, sep = "/"))  # reading temperature files
print("Loaded temperature rasters")
# TODO: Possibly needs change for other countries
ta <- ta / 10 # adjusting temperature values (originally multiplied by 10)
print("Loaded and adjusted temperature data")

pr <- list.files(pavg95)  # list of precipitation files
pr <- pr[grepl(".tif$", pr)]  # filtering only tif files
pr <- rast(paste(pavg95, pr, sep = "/"))  # reading precipitation files
# TODO: Possibly needs change for other countries
pr <- pr / 100 * 16 # adjusting precipitation values (originally multiplied by 100 and for 100m raster cells)
print("Loaded and adjusted precipitation data")

ph <- rast(ph)
ph <- aggregate(ph, fact = 4, fun = mean)
print("Loaded and aggregated pH data")

##---- Ecocrop model on each crop type
# counter variable
z <- 1
for (i in 1:nrow(crops)) {
  # Extract the crop name from the current row.
  a <- crops$crop[i]
  # Create a short version of the name containing 3 characters
  c_name <- substr(a, start = 1, stop = 3)
  # if short version allready exists in crops$mod then z is added to the name
  if (c_name %in% crops$mod_name == TRUE) {
    crops$mod_name[i] <- paste(c_name, z, sep = "_")
    z <- z + 1 }
  else {
    crops$mod_name[i] <- c_name #mod_name becomes a column with each model name
  }
  # Get parameters for the current crop and create a model
  pars <- ecocropPars(a)
  m_name <- ecocrop(pars) #model
  print(paste(
    "Iteration:", i,
    "Crop:", a,
    "Short name:", c_name,
    sep = " "
  ))
  control(m_name, get_max = TRUE)  # get the maximum value of the model
  print("Model created")
  # Predict using the model
  pred <- predict(m_name, tavg = ta, prec = pr, ph = ph, wopt = list(names = c_name))
  crs(pred) <- crs(params$proj$crs)
  # Export the prediction
  ex_name <- paste(c_name, ".tif", sep = "")
  terra::writeRaster(pred, file.path(out_dir, ex_name), overwrite = TRUE)
  # Print progress
  print(paste("model:", a, "done.", round((i * 100) / nrow(crops), 1), "%", sep = " "))
  gc()
}
