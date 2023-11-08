###################
#
# Aggregating species of interest of Medicinal resources NCP
#
# The script reads a list of species and their evaluations. 
# The species list is filtered to exclude aquatic species and further refined to include only species with 
# MED ES values. For each species in the refined list, 
# the script identifies its corresponding Species Distribution Model (SDM) map. 
# These SDM maps are then stacked together, and the mean and sum of the stacked maps are calculated
#
###################

# Load required libraries
library(terra)
library(data.table)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

start <- Sys.time()

# Output folder
out_dir <- params$MED$out_dir

# SDM map folder
sdm_curr_dir <- params$data$sdm_curr_dir

# Read species list and evaluations data
sp_list <- fread(params$MED$sp_list)
evals_list <- fread(params$MED$evals_list)

# Removing aquatic species
sp_list$terr <- NA

# For each species in the list, check its type (terrestrial or aquatic) from the evaluations list.
for (i in 1:nrow(sp_list)) {
  # Extract the binomial name of the species from the current row of the species list
  a <- sp_list$binom_name[i]
  # Find the row(s) in the evaluations list that matches the current species' binomial name
  b <- evals_list[evals_list$species == a,]

  # Check if there are any matching rows in the evaluations list for the current species
  if (nrow(b) == 0) {
    next
  } else {
    #If a match is found, extract the 'aqua_terr' value to the 'terr' column of the current row in the species list
    sp_list$terr[i] <- b$aqua_terr
  }
}

# Filter out aquatic species and further filter for species with specific ES values
sp_list <- sp_list[sp_list$terr != "aqua",]
list_sp <- sp_list[sp_list$ES == "MED" | sp_list$ES == "MEDPLR"]

print(paste("Number of species:", nrow(list_sp), sep = " "))

# Format species names and add filename column
list_sp <- as.data.frame(lapply(list_sp, function(y) gsub(" ", ".", y)))

# Get the list of SDM files with .tif extension
list_sdm <- list.files(sdm_curr_dir, pattern = "\\.tif$")

# For each species, find its corresponding SDM map filename
for (i in 1:nrow(list_sp)) {
  # For each species in the list, retrieve its name
  a <- list_sp$binom_name[i]
  # Search for the name in the list of SDM files.
  b <- grep(a, list_sdm)
  if (length(b) == 0) {
    next
  } else {
    # If found, assign the corresponding SDM filename to the 'filename' column of the species list
    list_sp$filename[i] <- list_sdm[b]
  }
}
# Construct the full path for each SDM map
list_sp$fullfilename <- paste(sdm_curr_dir, list_sp$filename, sep = "/")

# Remove rows with missing filenames and duplicates
list_sp <- na.omit(list_sp)
list_sp <- list_sp[-which(duplicated(list_sp$fullfilename)),]

# Stack all the SDM maps of the filtered species
stack_sdm_med <- rast(list_sp$fullfilename)

# Calculate the mean and sum of the stacked SDM maps
med_mean <- mean(stack_sdm_med)
print("Calculation of mean done")
med_sum <- sum(stack_sdm_med)
print("Calculation of sum done")

# Write rasters and list
writeRaster(med_mean, paste(out_dir, "med_mean.tif", sep = "/"), overwrite = TRUE)
writeRaster(med_sum, paste(out_dir, "med_sum.tif", sep = "/"), overwrite = TRUE)
# export a list of the species aggregated
fwrite(list_sp, paste(out_dir, "list_sp_mean.csv", sep = "/"))

print("Rasters and list written")

stop <-
  elapsed <- stop - start

# simple: print(elapsed)
# redable format:
print(paste(
  "Elapsed time: ",
  format(Sys.time() - start, units = "mins"),
  sep = " "
))