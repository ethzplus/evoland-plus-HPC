################
#
#The script reads a list of species and filters it to retain only those associated with pest control.
#It then maps each species to its corresponding SDM map filename.
#Then a raster stack is created and populated with SDM maps for each pest control species.
#Finally the sum and mean of those SDM maps is calculated
#
################

library(data.table)
library(terra)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

# Set working directory
setwd(params$PC$wd)

# SDM map folder
sdm_fold <- params$data$sdm_glo_dir

# output folder

outfold <- paste(wd, "output", sep = "/")

# Read species list for pest control
sp_list <- fread(params$PC$sp_list)

# Filter the species list to species associated with pest control (ES == "PC")
list_sp <- sp_list[sp_list$ES == "PC"]

# Format species names and add filename column
list_sp <- as.data.frame(lapply(list_sp, function(y) gsub(" ", ".", y)))
list_sp$filename <- NA

# Get a list of all SDM map filenames
list_sdm <- list.files(sdm_fold)

# For each species in the list, find its corresponding SDM map filename and store it in the 'filename' column
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

# Stack pest control species SDM maps
# create empty raster "stack"
stack_sdm_pc <- rast()
# Populate the raster stack with SDM maps for each species
for (i in 1:nrow(list_sp)) {
  # For each species in the list, retrieve the filename of its corresponding SDM map
  a <- list_sp$filename[i]
  #skip if filename is missing
  if (is.na(a) == TRUE) {
    next
  } else {
    # Construct the full path to the SDM map
    b <- rast(paste(sdm_fold, a, sep = "/"))
    # add SDM map to raster stack
    stack_sdm_pc <- c(stack_sdm_pc, b)
    print(a)
  }
}
# Remove rows with missing filenames from the species list
list_sp <- na.omit(list_sp)

# Calculate the sum and mean of the SDM maps in the raster stack
pc_1 <- sum(stack_sdm_pc)
pc_2 <- mean(stack_sdm_pc)

# Export to raster
#saveRDS(readAll(pc_1), "pest_control/FINAL/SUPPLY/18/sum_pc.rds")
#saveRDS(readAll(pc_2), "pest_control/FINAL/SUPPLY/18/mean_pc.rds")


# needs change, also writeRaster pc_1 if needed
writeRaster(pc_2, paste(outfold, "PC_S_CH_18.tif", sep = "/"))
# possibly wrong format list as .tif?
fwrite(list_sp, paste(outfold, "list_sp_mean.tif", sep = "/"))  #export a list of the species aggregated
