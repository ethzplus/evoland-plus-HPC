##################################?
#
# The script processes Species Distribution Models (SDM) maps for symbolic
# species.
#
# It starts by loading a list of species and filters it to retain only
# symbolic species. The script then matches each species with its
# corresponding SDM map file. After matching, it stacks the SDM maps of
# symbolic species and calculates both the mean and sum of the stacked maps.
#
##################################

library(terra)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

# Set working directory
wd <- params$ID$wd
setwd(wd)

# SDM map folder
sdm_reg_cov_dir <- params$data$sdm_reg_cov_dir
outdir_18 <- params$ID$outdir_18

# Loading sp list
list_sp <- fread(params$ID$list_sp, sep = ";")
list_sp <- list_sp[list_sp$ES == "SYM",] # Keep only symbolic species

# # Format the species names by replacing spaces with dots
list_sp <- as.data.frame(lapply(list_sp, function(y) gsub(" ", ".", y)))
# Initialize a new column to store filenames
list_sp$filename <- NA

# Get list of SDM files
list_sdm <- list.files(sdm_reg_cov_dir, pattern = "\\.tif$")

# Match species names with their corresponding SDM filenames
for (i in 1:nrow(list_sp)) {
  a <- list_sp$binom_name[i]
  b <- grep(a, list_sdm)

  # Check if the binomial name was found in the list of SDM files
  if (length(b) == 0) {
    next
  } else {
    # If the name was found, assign the corresponding filename from
    # 'list_sdm' to the 'filename' column of 'list_sp'
    list_sp$filename[i] <- list_sdm[b]
  }
}

# Remove rows with missing filenames
list_sp <- na.omit(list_sp)
# Construct the full path for each SDM file
list_sp$fullname <- paste(sdm_reg_cov_dir, list_sp$filename, sep = "/")

# Stacking symbolic species SDM maps
stack_sdm_pc <- rast(list_sp$fullname)
# Calculate the mean and sum of the stacked SDM maps
sym_1 <- mean(stack_sdm_pc)
sym_2 <- sum(stack_sdm_pc)

# Set CRS and write raster
crs(sym_1) <- params$proj$crs
# Depends on the output desired (mean or sum)
writeRaster(sym_1, paste(outdir_18, "ID_S_CH_18.tif", sep = "/"))
# export a list of the species aggregated
fwrite(list_sp, paste(outdir_18, "list_sp_mean.csv", sep = "/"))

