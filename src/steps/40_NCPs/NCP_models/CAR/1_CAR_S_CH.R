#####
#
# The aim of this script is to reclassify the Land-use map for the carbon
# storage NCP mapping, based on altitude (DEM) and production region of
# Switzerland.
# The outputs are one land-use map per Region/elevation.
#
# The scipt loads elevation (DEM) and land-use/land-cover (LULC) data,
# reclassifies the DEM into altitude classes, and intersects it with
# production regions.
# The script then reclassifies the LULC data into specific categories and
# creates LULC rasters for each region based on the intersection of elevation
# and production regions.
#
#####

# Load libraries
library(terra)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

# Digital elevation model
dem <- rast(params$data$dem)
# Production regions from CH
prodreg <- vect(params$data$prodreg)
# Land use
lulc18 <- rast(params$data$lulcc2018)

# Folder paths
setwd(params$CAR$wd)

# Creating scratch folders used for temporary storage of data
scratch <- paste(params$CAR$wd, "scratch", sep = "/")
dir.create(paste(params$CAR$wd, "scratch", sep = "/"))
dir.create(paste(scratch, "lulc_clip", sep = "/"))
dir.create(paste(scratch, "lulc_clip", "18", sep = "/"))
dir.create(paste(scratch, "Invest_models", sep = "/"))

# 1) Preprocessing data

# TODO: Possibly needs change for Peru
# Reclassify DEM in 3 altitude classes 
m1 <- c(
  0, 600, 1,
  600, 1200, 2,
  1200, Inf, 3
)
m11 <- matrix(m1, ncol = 3, byrow = TRUE)
dem_r <- classify(dem, m11)

# Convert raster to polygons
dem_p <- as.polygons(dem_r)

# Intersect elevation and production region
regelev <- intersect(dem_p, prodreg)
r <- regelev

# Create new column with region + elevation class in the attribute table of regelev
regelev$ProdregN_1 <- gsub("é", "e", regelev$ProdregN_1) # Remove accents
regelev$ProdregN_1 <- gsub(" ", "", regelev$ProdregN_1) # Remove tabs
regelev$regelev_n <- paste(regelev$ProdregN_1, regelev$DEM_mean_LV95, sep = "")

# Function to apply to each LULC map

# Reclassify the LULC map into 18 categories
# Always pairs of two, first number is original LULC category and second number
# the one it is casted into
# TODO: change to use the LULC-classes created by Ben
clip_and_reclassify <- function(lulc, year) {
  m <- c(
    0, 0, 10, 12, 11, 12, 12, 2, 13, 1, 14, 4, 15, 3, 16, 18, 17, 18, 7, 19, 16
  )

  #making the vector a matrix with two columns
  m1 <- matrix(m, ncol = 2, byrow = TRUE)
  lulc_r <- classify(lulc, m1)

  # Creating the rasters for each region
  list_reg_elev <- data.frame(unique(regelev$regelev_n))
  colnames(list_reg_elev) <- "regelev"

  # Iterates through each unique combination of production region and elevation class.
  for (i in 1:nrow(list_reg_elev)) {
    # Extracts the name of the current combination.
    name <- list_reg_elev$regelev[i]
    # Subsets, crops and masks the 'regelev' data for the current combination
    a <- regelev[regelev$regelev_n == name,]
    b <- crop(lulc_r, a)
    c <- mask(b, a)

    # Sets the NA (missing data) flag value to 255
    NAflag(c) <- 255

    # saving the data
    writeRaster(c, paste(scratch, "lulc_clip", year, paste(name, ".tif", sep = ""), sep = "/"), overwrite = TRUE, NAflag = 255)
    print(paste("raster", name, "created", i, "/", nrow(list_reg_elev), sep = " "))
    gc()
  }
}

# Applying the function
clip_and_reclassify(lulc18, "18")
