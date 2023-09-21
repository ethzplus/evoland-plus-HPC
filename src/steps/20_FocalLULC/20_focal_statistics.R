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
#' @author Carlson Büth, Benjamin Black
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

  # Split input raster into layer for each class
  seg_raster <- terra::segregate(raster)

  # Loop over layers calculating focal
  focal_layers <- terra::sapp(seg_raster, fun = function(lyr, ...) {
    terra::focal(lyr, focal_mat, fun = fun, ...)
  })

  # TODO: Name layers in Focal_layers according to LULC values or class names
  # (check require from Antoine A.)

  # This line would name according to raster values would need something else
  # to get class names
  names(focal_layers) <- unique(raster, na.rm = FALSE)

}

#' Convert map to predictor
#'
#' Loads the map from map_path, extracts the year, so focal_stats() can save
#' the predictor layer per class with the year in the name, in the format
#' "{base_path}_{year}_cl{class}_{radius}.rds"
#' The year is extracted from the map name by _(\d{4})(_|\.),
#' a four digit number enclosed by `_`, and `_` | `.` at the end.
#'
#' @param map_path Path to map, tif or rds.
#' Needs to include four digits for year
#' @param base_path Base path for predictor
#' @param ... Additional arguments for focal_stats
#'
#' @examples
#' map_to_predictor("SimulationMap_2020.tif", "folder/name", radius = 500)
#' # -> "folder/name_2025_cl1_500.rds", "folder/name_2025_cl2_500.rds", ...
#'
#' @docType methods
#' @importFrom terra rast
#' @importFrom stringr str_match
#'

map_to_predictor <- function(map_path, base_path, ...) {
  # Load map
  map <- terra::rast(map_path)
  # Match year in map name by _(\d{4})(_|\.) and extract first capture group
  year <- stringr::str_match(basename(map_path), "_(\\d{4})(_|\\.)")[, 2]
  base_path <- paste0(base_path, "_", year)

  # Calculate focal statistics
  focal_layers <- focal_stats(map, ...)

  rm(map)

  # focal_stats will now return a Spatraster with focal layers for each class
  # value loop over the layers saving an .rds for each and modifying the save
  # path.

  terra::sapp(focal_layers, fun = function(lyr, ...) {

    # modify save path: {base_path}_{year}_cl{class}_{radius}.rds
    layer_path <- paste0(
      base_path,
      "_cl", names(lyr),
      "_", list(...)$radius,
      ".rds"
    )

    # Save layer as rds
    saveRDS(lyr, layer_path)
  }) #close loop over layers

  # Free memory
  rm(focal_layers)
}

#' Convert all maps to predictors in folder
#'
#' @param folder Folder with simulation maps
#' @param save_folder Folder to save predictors
#' @param base_name Base string for naming convetion
#' @param combinations List of combinations of radius, windowType, fun, ...
#'
#' @examples
#' combinations <- list(
#'  list(radius = 100, window_type = "circle", fun = "mean"),
#'  list(radius = 500, window_type = "circle", fun = "mean", na.rm = TRUE),
#' )
#' maps_to_predictors("SimulationMaps", "Predictors", "tests", combinations)
#' # -> "Predictors/tests_2025_cl1_100.rds",
#' #    "Predictors/tests_2025_cl1_500.rds", ...
#'
#' @docType methods
#' @importFrom utils txtProgressBar setTxtProgressBar getTxtProgressBar
#'

folder_to_predictors <- function(folder, save_folder, base_name, combinations) {
  # Get all maps - tif or rds
  map_paths <- list.files(folder, full.names = TRUE)
  map_paths <- map_paths[grepl(".tif|.rds", map_paths)]
  # For each map
  progress <- utils::txtProgressBar(
    min = 0,
    max = length(map_paths),
    style = 3
  )
  cat(paste0("Calculating focal statistics for ", length(map_paths),
             " maps...\n"))
  for (map in map_paths) {
    # For each combination
    for (combination in combinations) {
      # Convert map to predictor
      do.call(
        map_to_predictor,
        c(
          list(map_path = map),
          list(base_path = file.path(save_folder, base_name)),
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
# Check if BaseName is set
if (is.null(config$BaseName) || config$BaseName == "") {
  stop("BaseName not set in config.json")
}

cat("Calculating focal statistics for LULC preparation\n")

# Set working directory from config$WorkDir
setwd(config$WorkDir)
cat("Working directory set to:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")
cat("Base name set to:", config$BaseName, "\n")
cat("Radii:", paste(config$RadiusList, collapse = ", "), "\n")


# Convert maps to predictors for each radius in config$RadiusList
folder_to_predictors(
  folder = config$InputDir,
  save_folder = config$OutputDir,
  base_name = config$BaseName,
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


# simulated_LULC_scenario_BAU_simID_v1_year_2020.tif
#           -> ch_lulc_agg11_future_pixel_2020_cl1_100.rds
#
# 'ch_lulc_geostat65_present_pixel_2013_2018_cl53_100.rds'
