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
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

#setting working directory
wd <- params$FF$ecocrop_rawdata_dir
setwd(wd)

#--local variables
# setting paths to local variables
sm_dir <- paste(wd, "results", sep = "/")
pavg95 <- paste(wd, "pavg_95/13_18/lv95", sep = "/") # precipitations
tavg95 <- paste(wd, "tavg_95/13_18/lv95", sep = "/") # temperature
ph <- paste(wd, "ch_edaphic_eiv_descombes_pixel_r.tif", sep = "/")      # pH
crops_data <- paste(wd, "crops.txt", sep = "/") #crops names

##---Loading variables

#List of food of interest

#read table containing crops list

crops <- t(fread(crops_data, sep = ";", header = FALSE))
rownames(crops) <- c(1:nrow(crops))
colnames(crops) <- "crop"
crops <- as.data.frame(crops)

#Environmental variables (Temp, prec, ph)

ta <- rast(paste(tavg95, list.files(tavg95), sep = "/"))
# TODO: Possibly needs change for other countries
ta <- ta / 10 # adjusting temperature values (originally multiplied by 10)

pr <- rast(paste(pavg95, list.files(pavg95), sep = "/"))
# TODO: Possibly needs change for other countries
pr <- pr / 100 * 16 # adjusting precipitation values (originally multiplied by 100 and for 100m raster cells)

ph <- rast(paste(ph, list.files(ph), sep = "/"))


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
  control(m_name, get_max = TRUE)
  # Predict using the model
  pred <- predict(m_name, tavg = ta, prec = pr, ph = ph, wopt = list(names = c_name))
  # Export the prediction
  ex_name <- paste(c_name, ".tif", sep = "")
  terra::writeRaster(pred, paste(sm_dir, ex_name, sep = "/"), overwrite = TRUE)
  # Print progress
  print(paste("model:", a, "done.", round((i * 100) / nrow(crops), 1), "%", sep = " "))
  gc()
}
