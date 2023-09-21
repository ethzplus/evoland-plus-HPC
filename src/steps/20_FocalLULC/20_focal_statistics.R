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
packs <- c("terra", "jsonlite", "future", "future.apply", "stringr")
invisible(lapply(packs, require, character.only = TRUE))


## Functions

#' Calculate focal statistics for a given raster and radius, and save as rds
#'
#' @param raster RasterLayer to calculate focal statistics for
#' @param radius Radius for window
#' @param base_path Base path for predictor
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
  raster, radius, base_path, window_type = "circle", fun = "mean", ...
) {
  # Check if radius is >= the resolution of the raster
  if (radius < min(terra::res(raster))) {
    stop(paste0(
      "Radius (", radius, ") is smaller than the resolution of the raster (",
      min(terra::res(raster)), ")."
    ))
  }
  # Create focal matrix
  focal_mat <- terra::focalMat(
    x = raster,
    d = radius,
    type = window_type
  )

  # Split input raster into layer for each class
  seg_raster <- terra::segregate(raster)

  # Loop over layers calculating focal
  sapply(as.list(seg_raster), function(lyr, ...) {
    focal_lyr <- terra::focal(
      x = lyr,
      w = focal_mat,
      fun = fun,
      ...
    )
    # modify save path: {base_path}_{year}_cl{class}_{radius}.rds
    layer_path <- paste0(
      base_path,
      "_cl", names(lyr),
      "_", radius,
      ".rds"
    )
    # Save layer as rds
    saveRDS(focal_lyr, layer_path)
    # Save to tif - terra::writeRaster(focal_lyr, paste0(layer_path, ".tif"))
    cat("Saved layer ", names(lyr), " to ", layer_path, "!\n")
  })

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

  # Skip if radius is smaller than resolution of map
  if (list(...)$radius < min(terra::res(map))) {
    cat(paste0(
      "Radius (", list(...)$radius, ") is smaller than the resolution of the ",
      "map (", min(terra::res(map)), "). Skipping map ", map_path, ".\n"
    ))
    return()
  }

  # Calculate focal statistics
  focal_stats(map, base_path, ...)
  rm(map)
  cat("Done calculating focal statistics for map ", map_path, "!\n")
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
#' @importFrom future future_sapply availableCores
#' @importFrom future.apply future_sapply
#'

folder_to_predictors <- function(
  folder, save_folder, base_name, combinations, parallel = FALSE
) {
  # Get all maps - tif or rds
  map_paths <- list.files(folder, full.names = TRUE)
  map_paths <- map_paths[grepl(".tif|.rds|.grd", map_paths)]
  # For each map
  if (parallel == FALSE) {
    cat(paste0("Calculating focal statistics for ", length(map_paths),
               " maps sequentially...\n"))
    progress <- utils::txtProgressBar(
      min = 0,
      max = length(map_paths) * length(combinations),
      style = 3
    )
    future::plan(future::sequential)
  } else {
    # number of workers from config - if not set, use all available cores
    n_workers <- ifelse(
      is.null(config$NWorkers) ||
        config$NWorkers == "" ||
        config$NWorkers == 0,
      future::availableCores(),
      config$NWorkers
    )
    cat(paste0("Calculating focal statistics for ", length(map_paths),
               " maps in parallel with ", n_workers, " workers...\n"))
    future::plan(future::multisession, workers = n_workers)
  }
  future.apply::future_sapply(map_paths, function(map) {
    cat(paste0("Calculating focal statistics for map ", map, "...\n"))
    # For each combination
    for (combination in combinations) {
      cat(paste0("Calculating focal statistics for combination ",
                 paste(combination, collapse = ", "), "...\n"))
      # Convert map to predictor
      do.call(
        map_to_predictor,
        c(
          list(map_path = map),
          list(base_path = file.path(save_folder, base_name)),
          combination
        )
      )
      if (parallel == FALSE) {
        # Update progress bar
        utils::setTxtProgressBar(
          progress,
          value = utils::getTxtProgressBar(progress) + 1
        )
      }
    }
  })
}


## Main

# Load config from src/config.json
config <- jsonlite::read_json("src/config.json", simplifyVector = TRUE)
config <- config$FocalLULCC # only FocalLULCC

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
  dir.create(config$OutputDir, recursive = TRUE)
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

t0 <- Sys.time()

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
        fun = config$FocalFunction
      )
    }
  ),
  parallel = config$Parallel
)

cat("Done calculating focal statistics for LULC preparation!\n")
t0 <- Sys.time() - t0
cat("Time elapsed: ", t0, attr(t0, "units"), "\n")
