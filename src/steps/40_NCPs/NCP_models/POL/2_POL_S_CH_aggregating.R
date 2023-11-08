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
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

#Local variables
setwd(params$POL$wd)

data_folder <- paste(getwd(), "/invest", sep = "")

# Retrieve a list of file paths from the data folder that match the pattern "supply"
files <- list.files(data_folder, pattern = "supply", full.names = T)

#output folder

out18 <- paste(getwd(), "18", sep = "/")


# Read and stack the rasters from the list of file paths into a single raster stack
p_list <- rast(files)

# Calculate the sum of the values across all rasters in the stack
pol_sum <- sum(p_list)

# Normalize the summed raster values to a range of 0 to 1
nx <- minmax(pol_sum)
rn <- (pol_sum - nx[1,]) / (nx[2,] - nx[1,])

writeRaster(rn, paste(out18, "POL_S_CH_18.tif", sep = "/"))

