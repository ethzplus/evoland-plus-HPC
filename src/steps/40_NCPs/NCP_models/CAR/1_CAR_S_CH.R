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

# Load the parameters into env by sourcing the ../load_params.R script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.dir <- dirname(sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]))
source(file.path(script.dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$data$dem)) { stop("params$data$dem is not set") }
if (is.null(params$data$prodreg)) { stop("params$data$prodreg is not set") }
if (is.null(params$run_params$NCP_RUN_SCENARIO_ID)) {
  stop("params$run_params$NCP_RUN_SCENARIO_ID is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$NCP_RUN_SCRATCH_DIR)) {
  stop("params$run_params$NCP_RUN_SCRATCH_DIR is not set")
}
if (is.null(params$data$lulc)) { stop("params$data$lulc is not set") }

# Digital elevation model
dem <- rast(params$data$dem)
# Production regions from CH
prodreg <- vect(params$data$prodreg)
# Land use
lulc <- rast(params$data$lulc)

# Creating scratch folder used for temporary storage of data
scratch <- file.path(params$run_params$NCP_RUN_SCRATCH_DIR,
                     params$run_params$NCP_RUN_SCENARIO_ID, "CAR")
out_dir <- file.path(scratch, "lulc_clip", params$run_params$NCP_RUN_YEAR)
if (dir.exists(out_dir))
  cat("Warning: Directory already exists, overwriting existing files.\n")
dir.create(
  out_dir,
  recursive = TRUE,  # creates parent directories if they don't exist
  showWarnings = FALSE
)

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
clip_and_reclassify <- function(lulc, out) {
  m <- c(
    0, 0,
    10, 12,
    11, 12,
    12, 2,
    13, 1,
    14, 4,
    15, 3,
    16, 18,
    17, 18,
    18, 7,
    19, 16
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
    filename <- file.path(out, paste(name, ".tif", sep = ""))
    writeRaster(c, filename, overwrite = TRUE, NAflag = 255)
    print(paste0(formatC(i, width = nchar(nrow(list_reg_elev)), flag = "0"),
                 "/", nrow(list_reg_elev), ": created ", name))
    gc()  # garbage collection
  }
}

# Applying the function
clip_and_reclassify(lulc, out_dir)