#' Focal statistics for LULC preparation
#'
#' Input: Folder name of simulation maps
#' for every simulation map:
#'     layerize grid to bands
#'     for each band:
#'         calculate focals
#'     save as rds
#' Output: N-SDM predictors
#'
#' @environment focal_lulc
#' @config config.json
#' @date 2023-09-20
#' @author Carlson Büth
#'
#' @docType script
#'

# Load libraries
packs <- c("terra", "jsonlite")
invisible(lapply(packs, require, character.only = TRUE))


## Functions

#' Calculate focal statistics for a given raster and radius
#'
#' @param raster RasterLayer to calculate focal statistics for
#' @param radius Radius for window
#' @param window_type Type of window, e.g. "circle", "rectangle"
#' @param fun Function to calculate focal statistics, function taking multiple
#' numbers and returning a numeric vector (one or more values),
#' e.g. "mean", "modal", "min", "max"
#' @param ... Additional arguments for fun, e.g. na.rm = TRUE
#' @return RasterLayer with focal statistics
#'
#' @examples
#' focal_stats(terra::rast(matrix(1:100, nrow = 10, ncol = 10)), 1)
#'
#' @docType methods
#' @importFrom terra focal focalMat
#'

focal_stats <- function(
  raster, radius, window_type = "circle", fun = "mean", ...
) {
  # Create focal matrix
  focal_mat <- terra::focalMat(
    x = raster,
    d = radius,
    type = window_type
  )
  # Calculate focal statistics
  terra::focal(raster, focal_mat, fun = fun, ...)
}

#' Convert map to predictor
#'
#' @param map_path Path to map, tif or rds
#' @param save_path Path to save predictor
#' @param ... Additional arguments for focal_stats
#'
#' @examples
#' map_to_predictor("SimulationMap_1.tif", "predictor_1.rds", radius = 500)
#'
#' @docType methods
#' @importFrom terra rast
#'

map_to_predictor <- function(map_path, save_path, ...) {
  # Load map
  map <- terra::rast(map_path)
  # Calculate focal statistics
  map <- focal_stats(map, ...)
  # Save as rds
  saveRDS(map, save_path)
  # Free memory
  rm(map)
}

#' Convert all maps to predictors in folder
#'
#' @param folder Folder with simulation maps
#' @param save_folder Folder to save predictors
#' @param combinations List of combinations of radius, windowType, fun, ...
#'
#' @examples
#' combinations <- list(
#'  list(radius = 100, window_type = "circle", fun = "mean"),
#'  list(radius = 500, window_type = "circle", fun = "mean", na.rm = TRUE),
#' )
#' maps_to_predictors("SimulationMaps", "Predictors", combinations)
#'
#' @docType methods
#'

maps_to_predictors <- function(folder, save_folder, combinations) {
  # Get all maps - tif or rds
  maps <- list.files(folder, full.names = TRUE)
  maps <- maps[grepl(".tif|.rds", maps)]
  # For each map
  progress <- utils::txtProgressBar(
    min = 0,
    max = length(maps),
    style = 3
  )
  cat(paste0("Calculating focal statistics for ", length(maps), " maps...\n"))
  for (map in maps) {
    # Get map name
    map_name <- basename(map)
    # For each combination
    for (combination in combinations) {
      # Get combination name: radius_windowType_fun_{mapName}.rds
      combination_name <- paste(
        combination$radius,
        combination$windowType,
        combination$fun,
        map_name,
        sep = "_"
      )
      # Convert map to predictor
      do.call(
        map_to_predictor,
        c(
          list(mapPath = map),
          list(savePath = file.path(save_folder, combination_name)),
          combination
        )
      )
    }
    # Update progress bar
    utils::setTxtProgressBar(
      progress,
      i = utils::getTxtProgressBar(progress) + 1
    )
  }
}


## Main

# Load config - only FocalLULCC
config <- jsonlite::read_json("../../config.json", simplifyVector = TRUE)
config <- config$FocalLULCC

# Check if InputDir is set
if (is.null(config$InputDir) || config$InputDir == "") {
  stop("InputDir not set in config.json")
}
# Check if OutputDir is set
if (is.null(config$OutputDir) || config$OutputDir == "") {
  stop("OutputDir not set in config.json")
}
# Check if OutputDir exists and create if not
if (!dir.exists(config$OutputDir)) {
  dir.create(config$OutputDir)
}

cat("Calculating focal statistics for LULC preparation\n")

# Set working directory from config$WorkDir
setwd(config$WorkDir)
cat("Working directory set to:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")
cat("Radii:", paste(config$RadiusList, collapse = ", "), "\n")


# Convert maps to predictors for each radius in config$RadiusList
maps_to_predictors(
  folder = config$InputDir,
  save_folder = config$OutputDir,
  combinations = lapply(
    config$RadiusList,
    function(radius) {
      list(
        radius = radius,
        windowType = config$WindowType,
        fun = config$FocalFunction,
      )
    }
  )
)

cat("Done calculating focal statistics for LULC preparation!\n")
