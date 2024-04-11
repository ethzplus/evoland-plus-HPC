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
#' @config $FUTURE_EI_CONFIG_FILE (yaml file)
#' @date 2023-09-20, 2023-12-08
#' @author Carlson Büth, Benjamin Black
#'
#' @docType script
#'

# Load libraries
packs <- c("terra", "raster", "yaml", "future", "future.apply", "stringr")
invisible(lapply(packs, require, character.only = TRUE))


## Functions

#' Calculate focal statistics for a given raster and radius, and save as rds
#'
#' @param raster RasterLayer to calculate focal statistics for
#' @param radius Radius for window
#' @param save_folder Folder to save predictors
#' @param base_name Base string for naming convetion
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
  raster, radius, save_folder, base_name, window_type = "circle", fun = "mean",
  ...
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
      na.rm = TRUE,
      na.policy = "omit",
      ...
    )
    layer_name <- paste0(base_name, "_reg", names(lyr), "_", radius)
    # Reformatting
    # Round to percent
    focal_lyr <- round(focal_lyr * 100)
    # Set storage mode to integer
    storage.mode(focal_lyr[]) <- "integer"
    # Name layer
    names(focal_lyr) <- layer_name
    # Convert back to raster
    focal_lyr <- raster::raster(focal_lyr)
    # convert file name to folder structure
    # a_b_c_d -> c(a, b, c, d) -> a/b/c (OS specific)
    layer_folder <- as.list(strsplit(layer_name, "_")[[1]])
    layer_folder <- do.call(
      file.path,
      layer_folder
    )
    save_folder <- file.path(save_folder, dirname(layer_folder))
    # Assure folder exists
    dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

    # Save layer as rds
    saveRDS(focal_lyr, file.path(save_folder, paste0(layer_name, ".rds")))
    # Save to tif - terra::writeRaster(focal_lyr, paste0(layer_path, ".tif"))
    cat("Saved layer ", names(lyr), " to ", layer_folder, "!\n")
  })

}

#' Convert map to predictor
#'
#' Loads the map from map_path, extracts the year, so focal_stats() can save
#' the predictor layer per class with the year in the name, in the format
#' "{save_folder}/{base_name}/{year}/{scenario_name}/reg{class}/
#' {base_name}_{year}_{scenario_name}_reg{class}_{radius}.rds"
#' The year is extracted from the map name by _(\d{4})(_|\.),
#' a four digit number enclosed by `_`, and `_` | `.` at the end.
#'
#' @param map_path Path to map, tif or rds.
#' Needs to include four digits for year
#' @param save_folder Folder to save predictors
#' @param base_name Base string for naming convetion
#' @param scenario_name Name of scenario
#' @param ... Additional arguments for focal_stats
#'
#' @examples
#' map_to_predictor("SimulationMap_2020.tif", "folder/name", radius = 500)
#' # -> "folder/name/2025/reg1/500.rds", "folder/name/2025/reg2/500.rds", ...
#'
#' @docType methods
#' @importFrom terra rast
#' @importFrom stringr str_match
#'

map_to_predictor <- function(
  map_path, save_folder, base_name, scenario_name, ...
) {
  # Load map
  map <- terra::rast(map_path)
  # Match year in map name by _(\d{4})(_|\.) and extract first capture group
  year <- stringr::str_match(basename(map_path), "_(\\d{4})(_|\\.)")[, 2]
  base_name <- paste0(base_name, "_", year, "_", scenario_name)

  # Skip if radius is smaller than resolution of map
  if (list(...)$radius < min(terra::res(map))) {
    cat(paste0(
      "Radius (", list(...)$radius, ") is smaller than the resolution of the ",
      "map (", min(terra::res(map)), "). Skipping map ", map_path, ".\n"
    ))
    return()
  }

  # Calculate focal statistics
  focal_stats(map, save_folder, base_name, ...)
  rm(map)
  cat("Done calculating focal statistics for map ", map_path, "!\n")
}

#' Convert all maps to predictors in folder
#'
#' @param folder Folder with simulation maps
#' @param save_folder Folder to save predictors
#' @param base_name Base string for naming convetion
#' @param scenario_name Name of scenario
#' @param combinations List of combinations of radius, windowType, fun, ...
#' @param parallel Whether to run in parallel
#'
#' @examples
#' combinations <- list(
#'  list(radius = 100, window_type = "circle", fun = "mean"),
#'  list(radius = 500, window_type = "circle", fun = "mean", na.rm = TRUE),
#' )
#' maps_to_predictors("SimulationMaps", "Predictors", "tests", combinations)
#' # -> "Predictors/tests/2025/reg1/100.rds",
#' #    "Predictors/tests/2025/reg1/500.rds", ...
#'
#' @docType methods
#' @importFrom utils txtProgressBar setTxtProgressBar getTxtProgressBar
#' @importFrom future future_sapply availableCores
#' @importFrom future.apply future_sapply
#'

folder_to_predictors <- function(
  folder, save_folder, base_name, scenario_name, combinations, parallel = FALSE
) {
  # Get all maps - tif or rds
  map_paths <- list.files(folder, full.names = TRUE)
  map_paths <- map_paths[grepl(".tif|.rds|.grd", map_paths)]
  # If no maps found, return
  if (length(map_paths) == 0) {
    cat(paste0("No maps found in ", folder, ". Skipping folder.\n"))
    return()
  }
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
  # Shuffle `map_paths` to distribute memory load
  future.apply::future_sapply(
    # nolint start: indentation_linter
    sample(map_paths),
    function(map) {
      cat(paste0("Calculating focal statistics for map ", map, "...\n"))
      # For each combination
      for (combination in combinations) {
        cat(paste0("Calculating focal statistics for combination ",
                   paste(combination, collapse = ", "), "...\n"))
        do.call( # Convert map to predictor
          map_to_predictor,
          c(
            list(map_path = map),
            list(save_folder = save_folder),
            list(base_name = base_name),
            list(scenario_name = scenario_name),
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
  # nolint end: indentation_linter
}

#' Convert all maps to predictors in subfolders
#'
#' Loops over all subfolders in folder and calls folder_to_predictors()
#' for each subfolder. Infers scenario name from subfolder name.
#'
#' If folder contains no subfolders, the folder itself will be processed.
#'
#' Folder structure:
#' folder
#' ├── Scenario1/v6/simulated_LULC_{...}_2020.tif
#' ├── Scenario1/v6/simulated_LULC_{...}_2025.tif
#' ├── Scenario1/v6/...
#' ├── Scenario2/v6/simulated_LULC_{...}_2020.tif
#' In this case, scenario name is "Scenario1" and "Scenario2". v6 is ignored.
#'
#' save_folder
#' ├── base/name/2020/Scenario1/reg1/base_name_2020_Scenario1_reg1_100.rds
#' ├── base/name/2020/Scenario2/reg1/base_name_2020_Scenario2_reg1_100.rds
#' ├── base/name/2025/Scenario1/reg1/base_name_2025_Scenario1_reg1_100.rds
#'
#' @param folder Folder with subfolders of simulation maps
#' @param save_folder Folder to save subfolder of predictors
#' @param base_name Base string for naming convetion
#' @param combinations List of combinations of radius, windowType, fun, ...
#' @param parallel Whether to run in parallel
#'
#' @docType methods

simulated_lulc_to_predictors <- function(
  folder, save_folder, base_name, combinations, parallel = FALSE
) {
  if (length(list.dirs(config$InputDir, recursive = FALSE)) == 0) {
    # Process only this folder
    cat(paste0("Calculating focal statistics only for ", folder, " ...\n"))
    folder_to_predictors(
      folder = folder,
      save_folder = save_folder,
      base_name = base_name,
      scenario_name = basename(folder),
      combinations = combinations,
      parallel = parallel
    )
    return()
  }
  # Get all subfolders
  subfolders <- list.dirs(folder, recursive = FALSE, full.names = TRUE)

  cat(paste0("Calculating focal statistics for ", length(subfolders),
             " subfolders in ", folder, " ...\n"))
  # For each subfolder
  for (subfolder in subfolders) {
    cat(paste0("Calculating focal statistics for subfolder ",
               subfolder, " ...\n"))
    # Convert maps to predictors
    folder_to_predictors(
      folder = subfolder,
      save_folder = save_folder,
      base_name = base_name,
      scenario_name = basename(subfolder),
      combinations = combinations,
      parallel = parallel
    )
  }
}


## Main

# Load config from $FUTURE_EI_CONFIG_FILE
config_file <- Sys.getenv("FUTURE_EI_CONFIG_FILE")
if (!file.exists(config_file)) {
  stop(paste0(
    "Config file FUTURE_EI_CONFIG_FILE (",
    config_file,
    ") does not exist."
  ))
}
config <- yaml.load_file(config_file)
bash_vars <- config$bash_variables # general bash variables
config <- config$FocalLULCC # only FocalLULCC

# if InputDir is not set, use bash variable LULCC_CH_OUTPUT_SIM_DIR
if (is.null(config$InputDir) || config$InputDir == "") {
  config$InputDir <- bash_vars$LULCC_CH_OUTPUT_SIM_DIR
}
# if OutputDir is not set, use bash variable FOCAL_OUTPUT_SIM_DIR
if (is.null(config$OutputDir) || config$OutputDir == "") {
  config$OutputDir <- bash_vars$FOCAL_OUTPUT_SIM_DIR
}

# Check if InputDir is now set
if (is.null(config$InputDir) || config$InputDir == "") {
  stop("InputDir nor LULCC_CH_OUTPUT_SIM_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if InputDir exists
if (!dir.exists(config$InputDir)) {
  stop("InputDir does not exist")
}
# Check if InputDir has subfolders
# - otherwise just this folder will be processed
if (length(list.dirs(config$InputDir, recursive = FALSE)) == 0) {
  # warn that only this folder will be processed
  cat("InputDir has no subfolders to process
 Only this folder will be processed\n")
}
# Check if OutputDir is set
if (is.null(config$OutputDir) || config$OutputDir == "") {
  stop("OutputDir nor FOCAL_OUTPUT_SIM_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if OutputDir exists and create if not
if (!dir.exists(config$OutputDir)) {
  dir.create(config$OutputDir, recursive = TRUE)
}
# Check if BaseName is set
if (is.null(config$BaseName) || config$BaseName == "") {
  stop("BaseName not set in FUTURE_EI_CONFIG_FILE")
}

cat("Calculating focal statistics for LULC preparation\n")

cat("Working directory set to:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")
cat("Base name set to:", config$BaseName, "\n")
cat("Radii:", paste(config$RadiusList, collapse = ", "), "\n")

t0 <- Sys.time()

# Convert maps to predictors for each radius in config$RadiusList
simulated_lulc_to_predictors(
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
