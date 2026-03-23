#' Calculate various summary metrics across the configurations
#'
#' Input: Folder of outputs from ES simulations:
#'    Normalise the ES outputs and calculate a series of different summary
#'    metrics from them, save results as rds files to be used in visualisations
#'    and subsequent analysis
#' Output: Summary measures of ES provision across configurations.
#'
#' @environment ES_summarisation
#' @config $FUTURE_EI_CONFIG_FILE (yaml file)
#' @date 2025-8-15
#' @author Benjamin Black
#'
#'
#' @docType script
#'

# # FOR TESTING ONLY
# # Create an environment variable for the config file
Sys.setenv(FUTURE_EI_CONFIG_FILE = "src/config.yml")
Sys.setenv(LULCC_M_SIM_CONTROL_TABLE = "Simulation_control.csv")

# Load libraries
packs <- c(
  "stringr",
  "terra",
  "future",
  "future.apply",
  "readxl",
  "data.table",
  "tidyr",
  "yaml",
  "dplyr",
  "viridis",
  "ggplot2",
  "tidyterra",
  "jsonlite",
  "magick",
  "grDevices",
  "classInt",
  "doParallel",
  "foreach",
  "sf",
  "raster",
  "biscale",
  "rlang" # Add this line
)

invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

### =========================================================================
### Functions
### =========================================================================

#' @describeIn util Ensure that a directory exists
ensure_dir <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  invisible(dir)
}

### Prepare a dataframe of information on the ES layers that are to be processed
#' Prepare a dataframe of information on the ES layers that are to be processed
#' Including the input data path for each config and time step,
#' a path to save the normalized tif file, a path to save a png file of the rescaled layer
#' a path to save a JSON file of the areas of land that fall within discrete classes of the ES provision values
#'
#' @param input_dir The directory where the input ES layers are stored
#' @param image_dir The directory where the output images will be saved
#' @param raster_dir The directory where the output rasters will be saved
#' @param perc_area_dir The directory where the classified area JSON files will be saved
#' @param base_dir The base directory where the output directories are located
#' @param Sim_ctrl_tbl_path The path to the simulation control table CSV file
#' @param ProjCH The projection CRS to use for the ES layers
#' @param ESs_to_summarise A character vector of ESs to summarise
#' @param ES_nesting A named list indicating whether each ES is nested (TRUE) or not (FALSE)
#' @param ES_file_names A named list of file names for each ES
#' @param map_masks A list of shapefiles to use as masks for the ES layers
#' @return A dataframe with information on the ES layers to be processed, including paths for input, output rasters, images, and classified areas
prepare_analysis_df <- function(
  input_dir = config$InputDir,
  image_dir = config$ImageDir,
  chg_image_dir = config$ChgImageDir,
  raster_dir = config$RasterDir,
  chg_raster_dir = config$ChgRasterDir,
  chg_raster_init = config$ChgRasterInit,
  sum_chg_raster_dir = config$SumChgRasterDir,
  perc_area_dir = config$PercAreaDir,
  area_chg_dir = config$AreaChgDir,
  sample_dir = "sample_data",
  base_dir = config$OutputDir,
  Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
  ProjCH = config$ProjCH,
  ESs_to_summarise = config$ESs_to_summarise,
  ES_nesting = config$ES_nesting,
  ES_file_names = config$ES_file_names,
  mask
) {
  # if the directories do not exist, create them
  ensure_dir(file.path(base_dir, image_dir))
  ensure_dir(file.path(base_dir, raster_dir))
  ensure_dir(file.path(base_dir, perc_area_dir))
  ensure_dir(file.path(base_dir, area_chg_dir))
  ensure_dir(file.path(base_dir, sample_dir))
  ensure_dir(file.path(base_dir, chg_image_dir))
  ensure_dir(file.path(base_dir, chg_raster_dir))
  ensure_dir(file.path(base_dir, chg_raster_init))
  ensure_dir(file.path(base_dir, sum_chg_raster_dir))

  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)

  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "EI",
    "EI_"
  )

  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "GREX",
    "GR_EX"
  )

  # convert all values in columns: Climate_scenario.string, Econ_scenario.string,
  # Pop_scenario.string and Scenario_ID.string to lower case
  Sim_ctrl_tbl$Climate_scenario.string <- tolower(
    Sim_ctrl_tbl$Climate_scenario.string
  )
  Sim_ctrl_tbl$Econ_scenario.string <- tolower(
    Sim_ctrl_tbl$Econ_scenario.string
  )
  Sim_ctrl_tbl$Pop_scenario.string <- tolower(Sim_ctrl_tbl$Pop_scenario.string)
  Sim_ctrl_tbl$Scenario_ID.string <- tolower(Sim_ctrl_tbl$Scenario_ID.string)

  # Get earliest scenario start date and latest end date
  Start_date <- min(Sim_ctrl_tbl$Scenario_start.real)
  End_date <- max(Sim_ctrl_tbl$Scenario_end.real)

  # Create seq of scenario time steps with Step_length.real
  Sim_time_steps <- seq(
    Start_date,
    End_date,
    by = Sim_ctrl_tbl$Step_length.real[1]
  )

  # Vector IDs of configurations to be analysed
  # Note Manually adjust this to analyse specific configurations
  Config_IDs <- unique(Sim_ctrl_tbl$Simulation_num.)

  # Loop over ESs_to_summarise and create vectors of file paths for their layers
  # according to Config_IDs, Sim_time_steps and nesting
  ES_paths <- lapply(ESs_to_summarise, function(ES) {
    # convert ES to lowercase
    ES_lower <- tolower(ES)

    # Loop over config_IDs returning paths as vector
    Config_paths <- lapply(Config_IDs, function(Config_ID) {
      # use Config_ID to subset Sim_ctrl_tbl to get scenario details
      Config_details <- Sim_ctrl_tbl[
        Sim_ctrl_tbl$Simulation_num. == Config_ID,
      ]

      # Create a vector of file paths for each time step
      Config_time_paths <- as.data.frame(rbindlist(lapply(
        Sim_time_steps,
        function(Time_step) {
          # If ES is not nested then append time step to file path
          if (ES_nesting[[ES]] == FALSE) {
            input_path <- file.path(
              input_dir,
              Config_ID,
              ES,
              paste0(ES_file_names[[ES]], "_", Time_step, ".tif")
            )
          } else if (ES_nesting[[ES]] == TRUE) {
            # If ES is nesting then append time step as a dir before file path
            input_path <- file.path(
              input_dir,
              Config_ID,
              ES,
              Time_step,
              paste0(ES_file_names[[ES]], ".tif")
            )
          }

          base_file_path <- paste0(
            ES_lower,
            "-",
            Time_step,
            "-",
            Config_details$Climate_scenario.string,
            "-",
            Config_details$Econ_scenario.string,
            "-",
            Config_details$Pop_scenario.string,
            "-",
            Config_details$Scenario_ID.string,
            "-",
            names(mask)
          )

          base_file_path_no_time <- paste0(
            ES_lower,
            "-",
            Config_details$Climate_scenario.string,
            "-",
            Config_details$Econ_scenario.string,
            "-",
            Config_details$Pop_scenario.string,
            "-",
            Config_details$Scenario_ID.string,
            "-",
            names(mask)
          )

          # tif path structure: NCP-Time_step-Config_details$Climate_scenario.string-Config_details$Econ_scenario.string-Config_details$Pop_scenario.string-Config_details$Scenario_ID.string-full.tif’
          tif_path <- file.path(
            base_dir,
            raster_dir,
            paste0(base_file_path, ".tif")
          )

          chg_tif_init_path <- file.path(
            base_dir,
            chg_raster_init,
            paste0(base_file_path, ".tif")
          )

          chg_tif_path <- file.path(
            base_dir,
            chg_raster_dir,
            paste0(base_file_path, ".tif")
          )

          sum_chg_tif_path <- file.path(
            base_dir,
            sum_chg_raster_dir,
            paste0(base_file_path_no_time, ".tif")
          )

          png_path <- file.path(
            base_dir,
            image_dir,
            paste0(base_file_path, ".png")
          )

          chg_png_path <- file.path(
            base_dir,
            chg_image_dir,
            paste0(base_file_path, ".png")
          )

          sample_path <- file.path(
            base_dir,
            sample_dir,
            paste0(base_file_path, ".rds")
          )

          perc_area_path <- file.path(
            base_dir,
            perc_area_dir,
            paste0(base_file_path, "-perc_area.json")
          )

          area_chg_path <- file.path(
            base_dir,
            area_chg_dir,
            paste0(base_file_path, "-perc_area_chg.json")
          )

          # combine all paths into a named list
          time_step_paths <- list(
            input_path,
            tif_path,
            chg_tif_init_path,
            chg_tif_path,
            sum_chg_tif_path,
            png_path,
            chg_png_path,
            sample_path,
            perc_area_path,
            area_chg_path
          )
          names(time_step_paths) <- c(
            "Path",
            "tif_path",
            "chg_tif_init_path",
            "chg_tif_path",
            "sum_chg_tif_path",
            "png_path",
            "chg_png_path",
            "sample_path",
            "perc_area_path",
            "area_chg_path"
          )
          return(time_step_paths)
        }
      ))) # close loop over time_steps

      #Add column for ES
      Config_time_paths$ES <- ES_lower

      # Add column for Config_ID
      Config_time_paths$Config_ID <- Config_ID

      # Add column for Time_step
      Config_time_paths$Time_step <- Sim_time_steps

      # loop over ES_time_paths and check which exist
      #Config_time_paths$Exists <- file.exists(Config_time_paths$Path)

      return(Config_time_paths)
    }) # close loop over Config_IDs

    # Bind list of dataframes into a single dataframe
    Config_df <- do.call(rbind, Config_paths)
  })

  #rbind list of dataframes into a single dataframe
  ES_path_df <- do.call(rbind, ES_paths)

  #check if all files exist
  # if (all(ES_path_df$Exists)) {
  #   message("All ES input files exist.")
  # } else {
  #   missing_layers <- ES_path_df[ES_path_df$Exists == FALSE,]
  #   stop("Some ES input files do not exist. Please check the paths.")
  # }

  return(ES_path_df)
}

#' calc_minmaxs:
#'
#' Calculate the minmax values for all configurations for each ES
#'
#' This function calculates the minimum and maximum values in the rasters
#' for all time points for all configurations for each ES  and saves
#' the results for each ES to an RDS file
#'
#'
#' @param parallel: logical: If TRUE, the function will run in parallel
#' using the number of workers specified in the config file. If FALSE, the
#' function will run sequentially
#' @param ES_layer_paths: data.frame: A df containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param ES_rescaling_dir: character: The directory to save the minmax
#' values for each ES
#' @param ESs_to_summarise: character: A vector of ESs to calculate minmax
#' values for
#' @param minmax_recalc: logical: If TRUE, the function will recalculate the
#' minmax values for each ES even if the file already exists. If FALSE, the
#' function will skip the calculation if the file already exists
#' @param report_NAs: logical: If TRUE, the function will check for NAs and NaNs
#' in the minmax values. If any are found, the function will return a vector of
#' the ES paths that have NAs or NaNs as the min or max values indicating
#' there may be a problem with the layers
#'
#'@return: vector: If report_NAs is TRUE, the function will return a vector of
#' the ES paths that have NAs or NaNs as the min or max values indicating
#' there may be a problem with the layers
#' @export
calc_minmaxs <- function(
  ES_layer_paths,
  map_path = "Path",
  save_path = ES_minmax_path,
  ES_rescaling_dir,
  report_NAs = TRUE,
  mask,
  ES
) {
  # Inner loop over paths for ES calculating the greatest minmax values for each
  ES_minmax <- future_sapply(
    ES_layer_paths[[map_path]],
    function(path) {
      # If mask is provided, apply it to the raster
      if (names(mask) != "full_extent") {
        # if the mask path contains shp extension, read it as a vector
        if (grepl("\\.shp$", unlist(mask))) {
          message(paste("Applying mask from shapefile:", mask))
          mask_layer <- terra::vect(x = unlist(mask))
        } else if (grepl("\\.tif$", mask)) {
          message(paste("Applying mask from raster file:", mask))
          mask_layer <- terra::rast(unlist(mask))
        } else {
          stop(paste("Unsupported mask file type for:", mask))
        }
      }

      # Calculate min and max for current ES
      # Use try catch in case the file is missing or corrupt
      File_minmax <- tryCatch(
        {
          # Load file
          raster <- rast(path)

          # If mask is provided, apply it to the raster
          if (names(mask) != "full_extent") {
            raster <- terra::crop(x = raster, y = mask_layer)

            raster <- terra::mask(
              x = raster,
              mask = mask_layer,
              updatevalue = NA
            )
          }

          # Get min max
          File_minmax <- minmax(raster, compute = TRUE)
        },
        error = function(e) {
          File_minmax <- NA
        }
      )

      cat("Finished calculating minmax for", path, "\n")
      return(File_minmax)
    },
    simplify = FALSE
  )

  # Save minmax values for ES
  saveRDS(ES_minmax, save_path)

  # if report_NAs is TRUE, check for NAs and NaNs in the minmax values
  if (report_NAs == TRUE) {
    # Summarise results across ESs in a single dataframe
    All_ES_minmaxs <- lapply(ESs_to_summarise, function(ES) {
      # Load the minmax values for current ES
      ES_minmax <- readRDS(file.path(
        Mask_rescaling_dir,
        paste0(tolower(ES), "_minmaxs.rds")
      ))

      #loop over list and add details to dataframe
      ES_minmax <- lapply(ES_minmax, function(x) {
        data.frame(Min = x[1], Max = x[2])
      })

      # bind list of dataframes into a single dataframe
      ES_minmax <- as.data.frame(rbindlist(ES_minmax, idcol = map_path))
    })
    names(All_ES_minmaxs) <- ESs_to_summarise

    # bind list of dataframes into a single dataframe
    All_ES_minmaxs <- as.data.frame(rbindlist(All_ES_minmaxs, idcol = "ES"))

    # Identify NAs
    NA_indices <- which(is.na(All_ES_minmaxs$Min))

    # Get paths of the NA value
    NA_records <- ES_layer_paths[NA_indices, map_path]

    # Check for NaN values
    NaN_indices <- which((is.nan(All_ES_minmaxs$Min)))

    # What is the path of the NaN value
    NaN_records <- ES_layer_paths[NaN_indices, map_path]
    return(c(NA_records, NaN_records))
  }
  return(ES_minmax)
}


#' calc_global_minmaxs:
#'
#' Calculate the global minmax values for each ES
#'
#' This function calculates global minimum and maximum values
#' (i.e. smallest minimum and greatest maximum value) for each ES across the
#' configurations*time points and saves the results to an RDS file
#'
#'
#' @param ESs_to_summarise Vector of ESs to summarise
#' @param ES_rescaling_dir Directory to load/save the results
#'
#' @export

calc_global_minmaxs <- function(
  ES_minmax = ES_minmax
) {
  #loop over list and add details to dataframe
  ES_minmax <- lapply(ES_minmax, function(x) {
    data.frame(Min = x[1], Max = x[2])
  })

  # bind list of dataframes into a single dataframe
  ES_minmax <- as.data.frame(rbindlist(ES_minmax))

  # Initialize min and max values
  ES_min <- Inf
  ES_max <- -Inf

  # Calculate global min and max for current ES
  for (i in 1:nrow(ES_minmax)) {
    ES_min <- min(ES_min, ES_minmax[i, "Min"], na.rm = TRUE)
    ES_max <- max(ES_max, ES_minmax[i, "Max"], na.rm = TRUE)
  }
  return(data.frame(Min = ES_min, Max = ES_max))
}


#' save_continuous_indexed_png
#' Save a continuous indexed PNG from a raster object with specified colors
#' called by normalise_layers function
#' @param raster_obj A RasterLayer or RasterStack object to be saved as PNG
#' @param output_path The file path where the PNG will be saved
#' @param low_color The color for the lowest value in the raster
#' @param high_color The color for the highest value in the raster
#' @param mid_color An optional color for the mid-point in the raster
#' @param width The width of the output PNG in pixels
#' @param height The height of the output PNG in pixels
#' @param resolution The resolution of the output PNG in DPI
#' @param units The units for width and height (default is "cm")
#' @param margins A numeric vector of length 4 specifying the margins (bottom, left, top, right) in inches
#' @param background The background color of the PNG (default is "transparent")
#' @param colorspace The colorspace for the indexed PNG (default is "sRGB")
#' @param max_colors The maximum number of colors to use in the indexed PNG (default is 256)
#' @param show_legend Whether to show the legend in the plot (default is FALSE)
#' @param axes Whether to show axes in the plot (default is FALSE)
#' @param box Whether to draw a box around the plot (default is FALSE)
#' @param cleanup_temp Whether to delete the temporary PNG file after saving (default is TRUE)
#' @param verbose Whether to print verbose messages during the process (default is FALSE)
save_continuous_indexed_png <- function(
  raster_obj,
  output_path,
  low_color,
  high_color,
  mid_color = NULL,
  low_value = 0, # New argument for low value
  high_value = 1, # New argument for high value
  mid_value = NULL, # New argument for mid value (used when mid_color is specified)
  width = 25,
  height = 20,
  resolution = 300,
  units = "cm",
  margins = c(0, 0, 0, 0),
  background = "transparent",
  colorspace = "sRGB",
  max_colors = 256,
  show_legend = FALSE,
  axes = FALSE,
  box = FALSE,
  cleanup_temp = TRUE,
  verbose = FALSE
) {
  # Load required libraries
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required but not installed.")
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package 'raster' is required but not installed.")
  }

  library(magick)
  library(raster)
  library(grDevices)

  # Validate color arguments
  if (missing(low_color) || missing(high_color)) {
    stop("You must provide at least 'low_color' and 'high_color'.")
  }

  # Validate value arguments
  if (low_value >= high_value) {
    stop("'low_value' must be less than 'high_value'.")
  }

  if (!is.null(mid_color)) {
    if (mid_value <= low_value || mid_value >= high_value) {
      stop("'mid_value' must be between 'low_value' and 'high_value'.")
    }
  }

  # Build color palette (with optional mid_color)
  if (!is.null(mid_color)) {
    if (verbose) {
      cat("Generating continuous palette with mid color...\n")
    }
    color_palette <- colorRampPalette(c(low_color, mid_color, high_color))(
      max_colors
    )
  } else {
    if (verbose) {
      cat("Generating continuous palette without mid color...\n")
    }
    color_palette <- colorRampPalette(c(low_color, high_color))(max_colors)
  }

  # Create temporary PNG
  temp_png <- tempfile(fileext = ".png")
  if (verbose) {
    cat("Creating temporary PNG...\n")
  }

  png(
    temp_png,
    width = width,
    height = height,
    res = resolution,
    bg = background,
    units = units
  )

  par(mar = margins)

  # Plot with explicit zlim to ensure consistent color mapping
  plot(
    raster_obj,
    col = color_palette,
    zlim = c(low_value, high_value), # This is the key addition
    legend = show_legend,
    axes = axes,
    box = box
  )

  dev.off()

  # Quantize to indexed PNG
  if (verbose) {
    cat("Converting to indexed PNG...\n")
  }
  img <- image_read(temp_png)

  img_indexed <- img %>%
    image_quantize(
      max = max_colors,
      colorspace = colorspace,
      dither = FALSE
    ) %>%
    image_strip()

  image_write(img_indexed, output_path, format = "png", depth = 8)

  if (verbose) {
    cat(paste("Saved indexed PNG to:", output_path, "\n"))
  }

  # Clean up
  if (cleanup_temp) {
    unlink(temp_png)
    if (verbose) cat("Cleaned up temporary files.\n")
  } else {
    if (verbose) cat(paste("Temporary file saved at:", temp_png, "\n"))
  }

  invisible(img_indexed)
}

#' normalise_layers:
#' Normalise the ES layers using the global minmax values
#' This function normalises the ES layers using the global minmax values
#' saving the layers as tif and png files using the provided paths in Es_layer_paths
#' called by normalise_summarize function
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param mask Character: The path to the mask to apply to the layers
#' @param ES_global_minmaxs Dataframe: containing the global min and max values for each ES
#' @return: NULL
normalise_layers <- function(
  ES_layer_paths,
  mask,
  ES_global_minmaxs,
  ProjCH = ProjCH,
  input_path_name = "Path",
  raster_path_name = "tif_path",
  image_path_name = "png_path",
  overwrite = TRUE,
  low_color = "#007CDC",
  high_color = "#FFEA2A",
  mid_color = NULL,
  low_value = 0,
  high_value = 1,
  mid_value = NULL
) {
  future_lapply(1:nrow(ES_layer_paths), function(i) {
    # Load the mask here because rast/vect are non-Future exportable objects
    if (names(mask) != "full_extent") {
      # if the mask path contains shp extension, read it as a vector
      if (grepl("\\.shp$", unlist(mask))) {
        message(paste("Applying mask from shapefile:", mask))
        mask_layer <- terra::vect(x = unlist(mask))
      } else if (grepl("\\.tif$", mask)) {
        message(paste("Applying mask from raster file:", mask))
        mask_layer <- terra::rast(unlist(mask))
      } else {
        stop(paste("Unsupported mask file type for:", mask))
      }
    }

    if (file.exists(ES_layer_paths[i, raster_path_name]) & overwrite == FALSE) {
      return(NULL)
    } else if (
      !file.exists(ES_layer_paths[i, raster_path_name]) | overwrite == TRUE
    ) {
      # If the input file doesn't exist skip
      if (!file.exists(ES_layer_paths[i, input_path_name])) {
        return(NULL)
      }

      # Load the raster
      raster <- rast(ES_layer_paths[i, input_path_name])

      # if the crs is not the same as ProjCH, reproject the raster
      if (crs(raster) != ProjCH) {
        raster <- terra::project(raster, ProjCH)
        cat("Reprojected raster to", ProjCH, "\n")
      }

      if (names(mask) != "full_extent") {
        # If mask is provided, apply it to the raster

        # crop raster to mask extent
        raster <- terra::crop(x = raster, y = mask_layer)

        raster <- terra::mask(x = raster, mask = mask_layer, updatevalue = NA)
      }

      ES_name <- ES_layer_paths[i, "ES"]
      ES_min <- ES_global_minmaxs[1, "Min"]
      ES_max <- ES_global_minmaxs[1, "Max"]

      # Normalise the raster
      raster_norm <- ((raster - ES_min) / (ES_max - ES_min)) *
        (high_value - low_value) +
        low_value

      #create dir
      dir <- dirname(ES_layer_paths[i, raster_path_name])
      ensure_dir(dir)

      # save as tif using the pre-prepared path
      writeRaster(
        raster_norm,
        file = ES_layer_paths[i, raster_path_name],
        overwrite = TRUE
      )

      cat("rescaled layer saved to", ES_layer_paths[i, raster_path_name], "\n")

      # produce png img and save
      save_continuous_indexed_png(
        raster_obj = raster_norm,
        output_path = ES_layer_paths[i, image_path_name],
        low_color = low_color,
        high_color = high_color,
        mid_color = mid_color,
        low_value = low_value, # New argument for low value
        high_value = high_value, # New argument for high value
        mid_value = mid_value, # New argument for mid value (used when mid_color is specified)
        width = 25,
        height = 20,
        resolution = 300,
        units = "cm",
        margins = c(0, 0, 0, 0),
        background = "transparent",
        colorspace = "sRGB",
        max_colors = 256,
        show_legend = FALSE,
        axes = FALSE,
        box = FALSE,
        cleanup_temp = TRUE,
        verbose = FALSE
      )
      cat(
        "rescaled layer image saved to",
        ES_layer_paths[i, image_path_name],
        "\n"
      )
    }
  })
}

#' normalise_layers:
#' Normalise the ES layers using the global minmax values
#' This function normalises the ES layers using the global minmax values
#' saving the layers as tif and png files using the provided paths in Es_layer_paths
#' called by normalise_summarize function
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param mask Character: The path to the mask to apply to the layers
#' @param ES_global_minmaxs Dataframe: containing the global min and max values for each ES
#' @return: NULL
normalise_layers_PATCH <- function(
  ES_layer_paths,
  mask,
  ES_global_minmaxs,
  ProjCH = ProjCH,
  input_path_name = "Path",
  raster_path_name = "tif_path",
  image_path_name = "png_path",
  overwrite = TRUE,
  low_color = "#007CDC",
  high_color = "#FFEA2A",
  mid_color = NULL,
  low_value = 0,
  high_value = 1,
  mid_value = NULL
) {
  future_lapply(1:nrow(ES_layer_paths), function(i) {
    # Load the mask here because rast/vect are non-Future exportable objects
    if (names(mask) != "full_extent") {
      # if the mask path contains shp extension, read it as a vector
      if (grepl("\\.shp$", unlist(mask))) {
        message(paste("Applying mask from shapefile:", mask))
        mask_layer <- terra::vect(x = unlist(mask))
      } else if (grepl("\\.tif$", mask)) {
        message(paste("Applying mask from raster file:", mask))
        mask_layer <- terra::rast(unlist(mask))
      } else {
        stop(paste("Unsupported mask file type for:", mask))
      }
    }

    if (file.exists(ES_layer_paths[i, raster_path_name]) & overwrite == FALSE) {
      return(NULL)
    } else if (
      !file.exists(ES_layer_paths[i, raster_path_name]) | overwrite == TRUE
    ) {
      # If the input file doesn't exist skip
      if (!file.exists(ES_layer_paths[i, input_path_name])) {
        return(NULL)
      }

      # # Load the raster
      # raster <- rast(ES_layer_paths[i, input_path_name])

      # # if the crs is not the same as ProjCH, reproject the raster
      # if (crs(raster) != ProjCH) {
      #   raster <- terra::project(raster, ProjCH)
      #   cat("Reprojected raster to", ProjCH, "\n")
      # }

      # if (names(mask) != "full_extent") {
      #   # If mask is provided, apply it to the raster

      #   # crop raster to mask extent
      #   raster <- terra::crop(x = raster, y = mask_layer)

      #   raster <- terra::mask(x = raster, mask = mask_layer, updatevalue = NA)
      # }

      # ES_name <- ES_layer_paths[i, "ES"]
      # ES_min <- ES_global_minmaxs[1, "Min"]
      # ES_max <- ES_global_minmaxs[1, "Max"]

      # # Normalise the raster
      # raster_norm <- ((raster - ES_min) / (ES_max - ES_min)) *
      #   (high_value - low_value) +
      #   low_value

      # #create dir
      # dir <- dirname(ES_layer_paths[i, raster_path_name])
      # ensure_dir(dir)

      # START BY RE-LOADING THE ALREADY NORMALISED RASTER TO AVOID RE-DOING THE NORMALISATION STEPS
      raster_norm <- rast(ES_layer_paths[i, raster_path_name])

      # mask IT TO THE NEW MASK IF NECESSARY
      if (names(mask) != "full_extent") {
        # If mask is provided, apply it to the raster

        # crop raster to mask extent
        raster_norm <- terra::crop(x = raster_norm, y = mask_layer)

        raster_norm <- terra::mask(
          x = raster_norm,
          mask = mask_layer,
          updatevalue = NA
        )
      }

      # save as tif using the pre-prepared path
      writeRaster(
        raster_norm,
        file = ES_layer_paths[i, raster_path_name],
        overwrite = TRUE
      )

      cat("rescaled layer saved to", ES_layer_paths[i, raster_path_name], "\n")

      # produce png img and save
      save_continuous_indexed_png(
        raster_obj = raster_norm,
        output_path = ES_layer_paths[i, image_path_name],
        low_color = low_color,
        high_color = high_color,
        mid_color = mid_color,
        low_value = low_value, # New argument for low value
        high_value = high_value, # New argument for high value
        mid_value = mid_value, # New argument for mid value (used when mid_color is specified)
        width = 25,
        height = 20,
        resolution = 300,
        units = "cm",
        margins = c(0, 0, 0, 0),
        background = "transparent",
        colorspace = "sRGB",
        max_colors = 256,
        show_legend = FALSE,
        axes = FALSE,
        box = FALSE,
        cleanup_temp = TRUE,
        verbose = FALSE
      )
      cat(
        "rescaled layer image saved to",
        ES_layer_paths[i, image_path_name],
        "\n"
      )
    }
  })
}

#' calc_summary_stats:
#' Calculate summary statistics for the rescaled ES layers
#' This function calculates the summary statistics for the rescaled ES layers
#' and returns the results in a dataframe
#' called by normalise_summarize function
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param metrics Vector of metrics to calculate must use built-in function arguments from terra::global
#' @return: A list of dataframes containing the summary statistics for each ES layer
calc_summary_stats <- function(
  ES_layer_paths,
  metrics = c("sum", "mean", "sd")
) {
  # Loop over sequence of ES layer paths
  ES_sum_stats <- future_lapply(1:nrow(ES_layer_paths), function(i) {
    # use trycatch to handle missing or corrupt files
    Layer_info_stat <- tryCatch(
      {
        # Load the already rescaled raster using the tif_path
        raster_norm <- rast(ES_layer_paths[i, "tif_path"])

        # Calculate the stats
        raster_stats <- global(raster_norm, metrics, na.rm = TRUE)

        # cbind raster_stats with the row of ES_layer_paths
        Layer_info_stat <- cbind(ES_layer_paths[i, ], raster_stats)

        # return the Layer_info_stat
        return(Layer_info_stat)
      },
      error = function(e) {
        # if an error occurs, return NA for all metrics
        raster_stats <- rep(NA, length(metrics))
        names(raster_stats) <- metrics

        #convert to dataframe
        raster_stats <- t(data.frame(raster_stats))

        # cbind raster_stats with the row of ES_layer_paths
        Layer_info_stat <- cbind(ES_layer_paths[i, ], raster_stats)
      }
    )
    return(Layer_info_stat)
  })
  return(ES_sum_stats)
}

#' tabulate_area_in_classes: Classify a raster against a single break vector and return bin percentages
#'
#' Extracted from Calculate_areas_in_classes so the same logic can be reused
#' for the whole-extent raster ("canton") and for each focus-area masked raster.
#' Returns a named list of bin_label -> percent suitable for embedding as one
#' entry inside a multi-area JSON object. JSON writing and plot saving are handled
#' by the caller so that all areas can be combined into a single output file.
#'
#' @param raster         A terra SpatRaster object to classify.
#' @param break_vector   Numeric vector of break boundaries for one break type.
#' @param context_label  String used in warning messages (e.g. "carbon – canton").
#' @return Named list of bin_label -> percent, or NULL if the raster has no data.
tabulate_area_in_classes <- function(raster, break_vector, context_label = "") {
  # ── Classify raster ──────────────────────────────────────────────────────
  rcl_matrix <- cbind(
    break_vector[-length(break_vector)],
    break_vector[-1],
    1:(length(break_vector) - 1)
  )
  rcl_matrix[1, 1] <- rcl_matrix[1, 1] - 0.001
  rcl_matrix[nrow(rcl_matrix), 2] <- rcl_matrix[nrow(rcl_matrix), 2] + 0.001

  r_classified <- terra::classify(
    raster,
    rcl = rcl_matrix,
    include.lowest = TRUE
  )

  # ── Frequency table ──────────────────────────────────────────────────────
  freq_tbl <- as.data.frame(freq(r_classified))
  freq_tbl$value[is.na(freq_tbl$value)] <- 0

  if (nrow(freq_tbl) == 0) {
    cat("Warning: No data for", context_label, "- skipping\n")
    return(NULL)
  }

  # ── Build complete bin table (fills zeros for missing bins) ──────────────
  n_expected_bins <- length(break_vector) - 1
  all_labels <- character(n_expected_bins)
  for (j in 1:n_expected_bins) {
    all_labels[j] <- paste0(
      round(break_vector[j], 2),
      " \u2013 ",
      round(break_vector[j + 1], 2)
    )
  }

  complete_freq_tbl <- data.frame(
    value = 1:n_expected_bins,
    count = 0L,
    stringsAsFactors = FALSE
  )

  for (k in 1:nrow(freq_tbl)) {
    bin_num <- freq_tbl$value[k]
    if (!is.na(bin_num) && bin_num >= 1 && bin_num <= n_expected_bins) {
      complete_freq_tbl$count[bin_num] <- freq_tbl$count[k]
    } else if (!is.na(bin_num)) {
      complete_freq_tbl <- rbind(
        complete_freq_tbl,
        data.frame(
          value = bin_num,
          count = freq_tbl$count[k],
          stringsAsFactors = FALSE
        )
      )
    }
  }

  total_cells <- sum(complete_freq_tbl$count)
  complete_freq_tbl$percent <- round(
    100 * complete_freq_tbl$count / total_cells,
    2
  )

  complete_freq_tbl$bin_label <- sapply(
    complete_freq_tbl$value,
    function(bin_num) {
      if (is.na(bin_num)) {
        return("NA")
      }
      if (bin_num < 1 || bin_num > length(all_labels)) {
        return(paste0("Bin_", bin_num, "_OutOfRange"))
      }
      all_labels[bin_num]
    }
  )

  # Return named list of bin_label -> percent (serialises as {"0 – 0.1": [84.42], ...})
  setNames(as.list(complete_freq_tbl$percent), complete_freq_tbl$bin_label)
}


# #' Calculate areas_in_classes:
# #' Calculate the areas of land that fall within discrete classes of the ES provision values
# #' This function calculates the areas of land that fall within discrete classes of the ES provision values
# #' and saves the results to a JSON file, it is capable of calculating using different break types including
# #' quantile, fisher and regular breaks that are calculated from samples taken acorss all the rescaled layers
# #' called by normalise_summarize function
# #' @param ES_layer_paths Dataframe: containing info for each ES layer
# #' i.e. config_ID, time step, raw layer path and path to save rescaled layer
# #' @param sample_size Integer: The number of samples to take from each layer for calculating global breaks
# #' @param break_types Character: The types of breaks to calculate, options are "quantile", "fisher" and "regular"
# #' @param save_samples Logical: Whether to save the samples to a file
# Calculate_areas_in_classes <- function(
#   ES_layer_paths,
#   sample_size = 50000,
#   save_samples,
#   break_types = c("quantile", "fisher", "regular"),
#   save_plots,
#   parallel,
#   n_workers = NULL,
#   input_map_path_name,
#   output_tbl_path_name,
#   regular_min = -1,
#   regular_max = 1,
#   regular_by = 0.25
# ) {
#   # break_types includes quantiles or fisher then it is necessary to take a
#   # sample from each raster layer in order to calculate global breaks under these methods
#   if (any(c("quantile", "fisher") %in% break_types)) {
#     if (parallel == TRUE) {
#       # Use provided n_workers or default to detectCores() - 1
#       if (is.null(n_workers) || n_workers == 0) {
#         n_workers <- detectCores() - 1
#       }

#       cl <- makeCluster(n_workers)
#       registerDoParallel(cl)

#       samples_list <- foreach(
#         i = 1:nrow(ES_layer_paths),
#         .combine = c,
#         .packages = "terra"
#       ) %dopar%
#         {
#           if (
#             save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])
#           ) {
#             samp <- readRDS(ES_layer_paths[i, "sample_path"])
#             return(samp)
#           } else if (file.exists(ES_layer_paths[i, "sample_path"]) == FALSE) {
#             raster_norm <- rast(ES_layer_paths[i, "tif_path"])
#             vals <- values(raster_norm, mat = FALSE)
#             vals <- vals[!is.na(vals)]
#             samp <- sample(
#               vals,
#               min(sample_size, length(vals)),
#               replace = FALSE
#             )
#             if (save_samples == TRUE) {
#               saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
#             }
#             return(samp)
#           }
#         }
#       stopCluster(cl)
#       all_samples <- samples_list
#     } else {
#       all_samples <- c()
#       for (i in 1:nrow(ES_layer_paths)) {
#         if (
#           save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])
#         ) {
#           samp <- readRDS(ES_layer_paths[i, "sample_path"])
#           all_samples <- c(all_samples, samp)
#           next
#         } else if (file.exists(ES_layer_paths[i, "sample_path"]) == FALSE) {
#           raster_norm <- rast(ES_layer_paths[i, input_map_path_name])
#           vals <- values(raster_norm, mat = FALSE)
#           vals <- vals[!is.na(vals)]
#           samp <- sample(vals, min(sample_size, length(vals)), replace = FALSE)
#           if (save_samples == TRUE) {
#             saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
#           }
#           all_samples <- c(all_samples, samp)
#         }
#       }
#     }
#   }

#   # default deciles for quantile/fisher
#   probs <- seq(0, 1, 0.1)

#   # empty breaks list
#   breaks <- list()
#   for (break_type in break_types) {
#     breaks[[break_type]] <- NA
#   }

#   for (break_type in break_types) {
#     if (break_type == "quantile") {
#       breaks$quantile <- unique(quantile(
#         all_samples,
#         probs = probs,
#         na.rm = TRUE
#       ))
#     } else if (break_type == "fisher") {
#       fisher <- classIntervals(all_samples, n = length(probs), style = "fisher")
#       breaks$fisher <- fisher$brks
#     } else if (break_type == "regular") {
#       breaks$regular <- seq(regular_min, regular_max, by = regular_by)
#     }
#   }

#   future_lapply(1:nrow(ES_layer_paths), function(i) {
#     raster <- rast(ES_layer_paths[i, input_map_path_name])

#     for (break_type in names(breaks)) {
#       current_breaks <- breaks[[break_type]]

#       if (break_type %in% c("quantile", "fisher")) {
#         rcl_matrix <- cbind(
#           current_breaks[-length(current_breaks)],
#           current_breaks[-1],
#           1:(length(current_breaks) - 1)
#         )
#         # Extend the range slightly to capture edge values
#         rcl_matrix[1, 1] <- rcl_matrix[1, 1] - 0.001 # Extend lower bound
#         rcl_matrix[nrow(rcl_matrix), 2] <- rcl_matrix[nrow(rcl_matrix), 2] +
#           0.001 # Extend upper bound
#         r_classified <- terra::classify(
#           raster,
#           rcl = rcl_matrix,
#           include.lowest = TRUE
#         )
#       } else if (break_type == "regular") {
#         rcl_matrix <- cbind(
#           current_breaks[-length(current_breaks)],
#           current_breaks[-1],
#           1:(length(current_breaks) - 1)
#         )
#         # Extend the range slightly to capture edge values
#         rcl_matrix[1, 1] <- rcl_matrix[1, 1] - 0.001 # Extend lower bound
#         rcl_matrix[nrow(rcl_matrix), 2] <- rcl_matrix[nrow(rcl_matrix), 2] +
#           0.001 # Extend upper bound
#         r_classified <- terra::classify(
#           raster,
#           rcl = rcl_matrix,
#           include.lowest = TRUE
#         )
#       }

#       freq_tbl <- freq(r_classified)
#       freq_tbl <- as.data.frame(freq_tbl)
#       freq_tbl$value[is.na(freq_tbl$value)] <- 0
#       if (nrow(freq_tbl) == 0) {
#         cat("Warning: No data for", break_type, "- skipping\n")
#         next
#       }

#       # Generate all expected bin numbers and labels
#       n_expected_bins <- length(current_breaks) - 1
#       all_bin_numbers <- 1:n_expected_bins
#       all_labels <- character(n_expected_bins)
#       for (j in 1:n_expected_bins) {
#         lower <- current_breaks[j]
#         upper <- current_breaks[j + 1]
#         all_labels[j] <- paste0(round(lower, 2), " – ", round(upper, 2))
#       }

#       # Create complete frequency table with all expected bins
#       complete_freq_tbl <- data.frame(
#         layer = rep(1, n_expected_bins), # Assuming single layer
#         value = all_bin_numbers,
#         count = 0, # Initialize all counts to 0
#         stringsAsFactors = FALSE
#       )

#       # Update counts for bins that actually have data
#       for (k in 1:nrow(freq_tbl)) {
#         bin_num <- freq_tbl$value[k]
#         if (!is.na(bin_num) && bin_num >= 1 && bin_num <= n_expected_bins) {
#           complete_freq_tbl$count[bin_num] <- freq_tbl$count[k]
#         } else if (!is.na(bin_num)) {
#           # Handle out-of-range bins by adding them to the complete table
#           out_of_range_row <- data.frame(
#             layer = 1,
#             value = bin_num,
#             count = freq_tbl$count[k],
#             stringsAsFactors = FALSE
#           )
#           complete_freq_tbl <- rbind(complete_freq_tbl, out_of_range_row)
#         }
#       }

#       # Calculate total cells for percentage calculation AFTER updating complete_freq_tbl
#       # This ensures out-of-range bins are included in the total
#       total_cells <- sum(complete_freq_tbl$count)

#       # Calculate percentages for the complete table
#       complete_freq_tbl$percent <- 100 * complete_freq_tbl$count / total_cells

#       # Create bin labels for all entries (including out-of-range)
#       complete_freq_tbl$bin_label <- sapply(
#         complete_freq_tbl$value,
#         function(bin_num) {
#           if (is.na(bin_num)) {
#             return("NA")
#           } else if (bin_num < 1 || bin_num > length(all_labels)) {
#             return(paste0("Bin_", bin_num, "_OutOfRange"))
#           } else {
#             return(all_labels[bin_num])
#           }
#         }
#       )

#       # Use the complete frequency table for further processing
#       freq_tbl <- complete_freq_tbl

#       freq_tbl$method <- break_type

#       # Round percentages to avoid floating point precision issues
#       freq_tbl$percent <- round(freq_tbl$percent, 2)

#       # Create unique keys by using bin_label as key instead of percent
#       # This ensures all bins are present and keys are unique
#       json_data <- toJSON(
#         setNames(as.list(freq_tbl$percent), freq_tbl$bin_label),
#         pretty = TRUE
#       )

#       # Write JSON to output path (no suffix since only using regular breaks)
#       write(json_data, file = ES_layer_paths[i, output_tbl_path_name])

#       if (save_plots == TRUE) {
#         chart_images_dir <- gsub(
#           "chart_data",
#           "chart_images",
#           dirname(ES_layer_paths[i, output_tbl_path_name])
#         )
#         ensure_dir(chart_images_dir)
#         chart_path <- file.path(
#           chart_images_dir,
#           gsub(
#             ".json",
#             ".png",
#             basename(ES_layer_paths[i, output_tbl_path_name])
#           )
#         )
#         bar_plot <- ggplot(freq_tbl, aes(x = bin_label, y = percent)) +
#           geom_bar(stat = "identity", fill = "steelblue") +
#           labs(
#             title = paste(
#               "Area in Classes for",
#               ES_layer_paths[i, "ES"],
#               "using",
#               break_type,
#               "breaks"
#             ),
#             x = "Class",
#             y = "Percentage of Area"
#           ) +
#           theme_minimal() +
#           theme(axis.text.x = element_text(angle = 45, hjust = 1))
#         ggsave(filename = chart_path, plot = bar_plot, width = 10, height = 6)
#       }
#     }
#   })
#   cat(
#     "Finished calculating breaks for all layers of current ES and saving classified area JSON files.\n"
#   )
# }

#' Calculate areas_in_classes:
#' Calculate the areas of land that fall within discrete classes of the ES provision values
#' This function calculates the areas of land that fall within discrete classes of the ES provision values
#' and saves the results to a JSON file, it is capable of calculating using different break types including
#' quantile, fisher and regular breaks that are calculated from samples taken across all the rescaled layers
#' called by normalise_summarize function
#' @param ES_layer_paths  Dataframe: containing info for each ES layer
#'   i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param sample_size     Integer: The number of samples to take from each layer for calculating global breaks
#' @param break_types     Character: The types of breaks to calculate, options are "quantile", "fisher" and "regular"
#' @param save_samples    Logical: Whether to save the samples to a file
#' @param focus_regions     Named list: Optional masks to apply to each layer before calculating area-in-classes.
#'   Each element is a file path to a raster that will be used as a mask.
#'   e.g. list("oberland_ost" = "Masks/focus_area_tifs/oberland_ost.tif")
#'   Results are written to a path derived from \code{output_tbl_path_name} by inserting
#'   the focus-area name before the file extension.
Calculate_areas_in_classes <- function(
  ES_layer_paths,
  sample_size = 50000,
  save_samples,
  break_types = c("quantile", "fisher", "regular"),
  save_plots,
  parallel,
  n_workers = NULL,
  input_map_path_name,
  output_tbl_path_name,
  regular_min = -1,
  regular_max = 1,
  regular_by = 0.25,
  focus_regions = NULL # <-- NEW: named list(area_name = mask_tif_path)
) {
  # ── Global break calculation (quantile / fisher need a pooled sample) ─────
  if (any(c("quantile", "fisher") %in% break_types)) {
    if (parallel == TRUE) {
      if (is.null(n_workers) || n_workers == 0) {
        n_workers <- detectCores() - 1
      }

      cl <- makeCluster(n_workers)
      registerDoParallel(cl)

      samples_list <- foreach(
        i = 1:nrow(ES_layer_paths),
        .combine = c,
        .packages = "terra"
      ) %dopar%
        {
          if (
            save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])
          ) {
            return(readRDS(ES_layer_paths[i, "sample_path"]))
          } else if (!file.exists(ES_layer_paths[i, "sample_path"])) {
            raster_norm <- rast(ES_layer_paths[i, "tif_path"])
            vals <- values(raster_norm, mat = FALSE)
            vals <- vals[!is.na(vals)]
            samp <- sample(
              vals,
              min(sample_size, length(vals)),
              replace = FALSE
            )
            if (save_samples == TRUE) {
              saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
            }
            return(samp)
          }
        }
      stopCluster(cl)
      all_samples <- samples_list
    } else {
      all_samples <- c()
      for (i in 1:nrow(ES_layer_paths)) {
        if (
          save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])
        ) {
          all_samples <- c(
            all_samples,
            readRDS(ES_layer_paths[i, "sample_path"])
          )
          next
        } else if (!file.exists(ES_layer_paths[i, "sample_path"])) {
          raster_norm <- rast(ES_layer_paths[i, input_map_path_name])
          vals <- values(raster_norm, mat = FALSE)
          vals <- vals[!is.na(vals)]
          samp <- sample(vals, min(sample_size, length(vals)), replace = FALSE)
          if (save_samples == TRUE) {
            saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
          }
          all_samples <- c(all_samples, samp)
        }
      }
    }
  }

  # ── Build breaks list ─────────────────────────────────────────────────────
  probs <- seq(0, 1, 0.1)
  breaks <- setNames(vector("list", length(break_types)), break_types)

  for (break_type in break_types) {
    if (break_type == "quantile") {
      breaks$quantile <- unique(quantile(
        all_samples,
        probs = probs,
        na.rm = TRUE
      ))
    } else if (break_type == "fisher") {
      fisher <- classIntervals(all_samples, n = length(probs), style = "fisher")
      breaks$fisher <- fisher$brks
    } else if (break_type == "regular") {
      breaks$regular <- seq(regular_min, regular_max, by = regular_by)
    }
  }

  # ── Per-layer processing ──────────────────────────────────────────────────
  future_lapply(1:nrow(ES_layer_paths), function(i) {
    raster <- rast(ES_layer_paths[i, input_map_path_name])
    base_output <- ES_layer_paths[i, output_tbl_path_name]
    es_label <- ES_layer_paths[i, "ES"]

    # Pre-load masked rasters once per layer (avoids re-reading inside break_type loop)
    masked_rasters <- list()
    if (!is.null(focus_regions) && length(focus_regions) > 0) {
      for (area_name in names(focus_regions)) {
        mask_path <- focus_regions[[area_name]]
        if (!file.exists(mask_path)) {
          cat(
            "Warning: mask file not found for focus area '",
            area_name,
            "': ",
            mask_path,
            " – skipping.\n",
            sep = ""
          )
          next
        }
        masked_rasters[[area_name]] <- terra::mask(raster, rast(mask_path))
      }
    }

    for (break_type in names(breaks)) {
      break_vector <- breaks[[break_type]]

      # ── Accumulate results: whole input raster first, then each focus area ─────────────
      # Output JSON structure:
      # {
      #   "canton":       { "0 – 0.1": [84.42], ... },
      #   "oberland_ost": { "0 – 0.1": [72.3],  ... },
      #   ...
      # }
      area_results <- list()

      area_results[["canton"]] <- tabulate_area_in_classes(
        raster = raster,
        break_vector = break_vector,
        context_label = paste(es_label, break_type, "canton")
      )

      for (area_name in names(masked_rasters)) {
        area_results[[area_name]] <- tabulate_area_in_classes(
          raster = masked_rasters[[area_name]],
          break_vector = break_vector,
          context_label = paste(es_label, break_type, area_name)
        )
      }

      # Drop any areas that returned NULL (no data / skipped)
      area_results <- Filter(Negate(is.null), area_results)

      browser()
      # ── Write combined JSON ─────────────────────────────────────────────────
      write(toJSON(area_results, pretty = TRUE), file = base_output)

      # ── Optional bar charts (one per area) ─────────────────────────────────
      if (save_plots) {
        chart_images_dir <- gsub(
          "chart_data",
          "chart_images",
          dirname(base_output)
        )
        ensure_dir(chart_images_dir)

        for (area_name in names(area_results)) {
          # Build a data frame from the bin list for ggplot
          bin_data <- data.frame(
            bin_label = names(area_results[[area_name]]),
            percent = unlist(area_results[[area_name]]),
            stringsAsFactors = FALSE
          )
          chart_path <- file.path(
            chart_images_dir,
            gsub(".json", paste0("_", area_name, ".png"), basename(base_output))
          )
          bar_plot <- ggplot(bin_data, aes(x = bin_label, y = percent)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            labs(
              title = paste(
                "Area in Classes –",
                es_label,
                "/",
                area_name,
                "(",
                break_type,
                "breaks)"
              ),
              x = "Class",
              y = "Percentage of Area"
            ) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          ggsave(filename = chart_path, plot = bar_plot, width = 10, height = 6)
        }
      }
    }
  })

  cat(
    "Finished calculating breaks for all layers of current ES and saving classified area JSON files.\n"
  )
}

#' Produce_change_map:
#' Produce change maps for each configuration
#' This function produces change maps for each configuration
#' and saves the results to the specified directories.
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param output_path_name Character: The name of the column in ES_layer_paths containing the output paths for the change maps
#' @param ES_chg_minmax_path Character: The path to save the minmax values for the change maps
#' @return: A list of minmax values for each change map
#' @export
Produce_change_map <- function(
  ES_layer_paths,
  diff_vs_init_path_name = "chg_tif_init_path", # column name for ti - t1 rasters
  sum_chg_path_name = "sum_chg_tif_path", # column name for summed-change raster
  ES_chg_minmax_path
) {
  # loop over each unique configuration ID for this scenario
  All_minmaxs <- future_lapply(
    unique(ES_layer_paths$Config_ID),
    function(Config_ID) {
      # subset the data to this Config_ID
      ES_config_paths <- ES_layer_paths[ES_layer_paths$Config_ID == Config_ID, ]

      cat(paste(
        "Calculating change in ES provision for Config_ID:",
        Config_ID,
        "\n"
      ))

      # order time steps and build ordered raster stack/SpatRaster
      time_order <- sort(unique(ES_config_paths$Time_step))
      paths_all <- ES_config_paths$tif_path[match(
        time_order,
        ES_config_paths$Time_step
      )]

      # if fewer than 2 time steps, nothing to do
      if (length(paths_all) < 2) {
        cat("Not enough time steps for Config_ID:", Config_ID, "- skipping\n")
        return(list())
      }

      # load all normalized rasters in chronological order
      ES_rast_all <- rast(paths_all)
      names(ES_rast_all) <- as.character(time_order)

      # first time step and subsequent comparison steps
      first_time_step <- time_order[1]
      comp_time_steps <- time_order[time_order != first_time_step]

      ES_rast_first <- ES_rast_all[[1]]
      ES_rast_comp_steps <- ES_rast_all[[2:nlyr(ES_rast_all)]]

      # calculate per-time-step differences vs first (ti - t1)
      ES_diff_from_first <- ES_rast_comp_steps - ES_rast_first
      names(ES_diff_from_first) <- as.character(comp_time_steps)

      # create an empty list for capturing the minmax values
      Time_step_minmaxs <- list()

      # loop over the time steps saving the image and calculating the area with different values of change in provision.
      for (time_step in as.character(comp_time_steps)) {
        cat(paste(
          "Saving image of change from",
          first_time_step,
          "to",
          time_step,
          "\n"
        ))

        Time_step_path <- ES_config_paths[
          as.character(ES_config_paths$Time_step) == as.character(time_step),
          diff_vs_init_path_name
        ]

        # ensure directory exists
        ensure_dir(dirname(Time_step_path[1]))

        # save as tif using the pre-prepared path
        writeRaster(
          ES_diff_from_first[[time_step]],
          file = Time_step_path[1],
          overwrite = TRUE
        )

        mm <- minmax(ES_diff_from_first[[time_step]], compute = TRUE)

        # calculate the min max values and add to the list
        Time_step_minmaxs[[time_step]] <- mm
      } # close for loop over time steps

      # --- Add cumulative sum of successive differences (same approach as calc_change_summaries) ---
      Diffs_over_time <- diff(ES_rast_all) # successive differences: t2-t1, t3-t2, ...
      Change_sum <- sum(Diffs_over_time) # cumulative sum across successive diffs

      # single save path (same for all time steps) taken from the specified column
      sum_chg_paths <- unique(ES_config_paths[[sum_chg_path_name]])
      if (length(sum_chg_paths) >= 1 && nzchar(sum_chg_paths[1])) {
        sum_chg_path <- sum_chg_paths[1]

        # ensure dir exists
        ensure_dir(dirname(sum_chg_path))

        # save summed change raster
        writeRaster(Change_sum, file = sum_chg_path, overwrite = TRUE)

        # compute minmax for summed change and include in returned list
        #mm_sum <- minmax(Change_sum, compute = TRUE)
        #Time_step_minmaxs[["sum_change"]] <- mm_sum

        cat("Saved sum-of-change raster to", sum_chg_path, "\n")
      } else {
        warning(
          "No sum_chg path found (column '",
          sum_chg_path_name,
          "') for Config_ID: ",
          Config_ID
        )
      }
      # --- END cumulative sum logic ---

      cat("finished all time steps for Config_ID:", Config_ID, "\n")
      return(Time_step_minmaxs)
    }
  ) # close for loop over unique Config_IDs

  All_minmaxs <- unlist(All_minmaxs, recursive = FALSE)

  # save the result
  saveRDS(All_minmaxs, ES_chg_minmax_path)
  return(All_minmaxs)
}

#' Undesirable_deviation
#'
#' Helper function used by Summarise_spatial_outputs to calculates the
#' undesirable deviation metric (Kwakkel et al. 2016) for a given set of values.
#' The undesirable deviation metric is a measure of robustness that incorporates
#' both a measure of performance and variation variance by calculating the
#' sum of the deviation (regret) away from the median performance value for all
#' instances below the median.
##'
#' @param f A numeric vector of values for which to calculate the undesirable deviation metric
#'
#' @references Kwakkel, Jan H., Sibel Eker, and Erik Pruyt. 2016.
#' ‘How Robust Is a Robust Policy? Comparing Alternative Robustness Metrics for
#'  Robust Decision-Making’. In Robustness Analysis in Decision Aiding,
#'  Optimization, and Analytics, edited by Michael Doumpos,
#'  Constantin Zopounidis, and Evangelos Grigoroudis, 221–37.
#'  Cham: Springer International Publishing.
#'  https://doi.org/10.1007/978-3-319-33121-8_10.
Undesirable_deviation <- function(f) {
  # calculating regret from median
  median <- median(f)

  # subtract the median from the value
  regret <- f - median

  # Now work with worst half of values

  # sort regret values
  regret <- sort(regret)

  # get the number of values
  n <- length(regret)

  # get the number of half of the values rounded up to nearest integer
  n_half <- round(n / 2 + 0.51)

  # get the worst half of the values
  worst_regret <- regret[1:n_half]

  # get the sum of the worst half of the values
  sum_worst_regret <- sum(worst_regret)

  return(sum_worst_regret)
}

#' summarize_cumulative_ES_change_maps:
#' Summarize cumulative ES change maps
#' This function summarizes cumulative ES change maps by calculating
#' spatial summary metrics such as mean, standard deviation,
#' mean-variance, and undesirable deviation across all ES change maps.
#' The summarized maps are saved to the specified directory.
#' @param sum_chg_rast_dir Directory containing the cumulative ES change maps
#' @param robustness_dir Directory to save the summarized maps
#' @param ES_groups A list of ES groups to process, each group is a vector of ES names
#' @param spatial_summary_metrics Vector of spatial summary metrics to calculate
#' @return: NULL
#' @export
summarize_cumulative_ES_change_maps <- function(
  sum_chg_rast_dir = config$SumChgRasterDir,
  robustness_dir = config$RobustnessDir,
  ES_groups = list(
    "ES" = c('REC', 'CAR', 'NDR', 'POL', 'SDR', 'WY', 'FF', 'HAB')
  ), # Changed from ESs_to_visualise
  spatial_summary_metrics = c(
    "mean",
    #"stdev",
    #"mean-variance",
    "undesirable-deviation"
  )
) {
  ensure_dir(robustness_dir)
  # Loop over each ES group
  for (group_name in names(ES_groups)) {
    ESs_to_visualise <- ES_groups[[group_name]]

    cat(
      "Processing ES group:",
      group_name,
      "with ES:",
      paste(ESs_to_visualise, collapse = ", "),
      "\n"
    )

    # List all tif raster paths
    rast_paths <- list.files(
      sum_chg_rast_dir,
      pattern = ".tif",
      full.names = TRUE
    )

    # subset to only those that contain the ES to visualise
    rast_paths <- rast_paths[grep(
      paste(ESs_to_visualise, collapse = "|"),
      rast_paths,
      ignore.case = TRUE
    )]

    # Skip if no rasters found for this group
    if (length(rast_paths) == 0) {
      cat("No rasters found for ES group:", group_name, "- skipping\n")
      next
    }

    # stack all of the rasters
    rast_stk <- rast(rast_paths)

    if (any(c("mean", "mean-variance") %in% spatial_summary_metrics)) {
      # calculate the mean across the raster layers
      rast_mean <- mean(rast_stk)
      mean_min <- minmax(rast_mean)[1]
      mean_max <- minmax(rast_mean)[2]

      # rescale values according to min max
      rast_mean_rescaled <- (rast_mean - mean_min) / (mean_max - mean_min)

      # Save the mean result with group name
      writeRaster(
        rast_mean,
        file.path(
          robustness_dir,
          paste0("Mean_sum_of_change_", group_name, ".tif")
        ),
        overwrite = TRUE
      )
      writeRaster(
        rast_mean_rescaled,
        file.path(
          robustness_dir,
          paste0("Mean_sum_of_change_", group_name, "_rescaled.tif")
        ),
        overwrite = TRUE
      )
    }

    if (any(c("stdev", "mean-variance") %in% spatial_summary_metrics)) {
      # calculate the standard deviation
      rast_stdev <- stdev(rast_stk)
      stdev_min <- minmax(rast_stdev)[1]
      stdev_max <- minmax(rast_stdev)[2]

      # rescale values according to min max
      rast_stdev_rescaled <- (rast_stdev - stdev_min) / (stdev_max - stdev_min)

      # Save the stdev result with group name
      writeRaster(
        rast_stdev,
        file.path(
          robustness_dir,
          paste0("Stdev_sum_of_change_", group_name, ".tif")
        ),
        overwrite = TRUE
      )
      writeRaster(
        rast_stdev_rescaled,
        file.path(
          robustness_dir,
          paste0("Stdev_sum_of_change_", group_name, "_rescaled.tif")
        ),
        overwrite = TRUE
      )
    }

    if ("mean-variance" %in% spatial_summary_metrics) {
      # Calculate the mean variance,the +1 is included to handle situations where
      # either the mean or standard deviation is close to zero.
      rast_mean_var <- (rast_mean + 1) / (rast_stdev + 1)

      # rescale values according to min max
      mean_var_min <- minmax(rast_mean_var)[1]
      mean_var_max <- minmax(rast_mean_var)[2]

      rast_mean_var_rescaled <- (rast_mean_var - mean_var_min) /
        (mean_var_max - mean_var_min)

      # Save the mean variance result with group name
      writeRaster(
        rast_mean_var,
        file.path(
          robustness_dir,
          paste0("Mean_variance_sum_of_change_", group_name, ".tif")
        ),
        overwrite = TRUE
      )
      writeRaster(
        rast_mean_var_rescaled,
        file.path(
          robustness_dir,
          paste0("Mean_variance_sum_of_change_", group_name, "_rescaled.tif")
        ),
        overwrite = TRUE
      )
    }

    # Do a clean up of already computed layers to free up memory
    suppressWarnings(
      rm(
        rast_mean,
        rast_mean_rescaled,
        rast_stdev,
        rast_stdev_rescaled,
        rast_mean_var,
        rast_mean_var_rescaled
      )
    )

    if ("undesirable-deviation" %in% spatial_summary_metrics) {
      message(
        "Calculating undesirable deviation spatial summary metric for group: ",
        group_name
      )

      # Calculate the undesirable deviation applying the function to all cells in the Spatraster
      terraOptions(nthreads = detectCores())
      rast_undesirable_dev <- app(rast_stk, fun = Undesirable_deviation)

      # rescale values according to min max
      undesirable_dev_min <- minmax(rast_undesirable_dev)[1]
      undesirable_dev_max <- minmax(rast_undesirable_dev)[2]

      rast_undesirable_dev_rescaled <- (rast_undesirable_dev -
        undesirable_dev_min) /
        (undesirable_dev_max - undesirable_dev_min)

      # Save the undesirable deviation result with group name
      writeRaster(
        rast_undesirable_dev,
        file.path(
          robustness_dir,
          paste0("Undesirable_deviation_sum_of_change_", group_name, ".tif")
        ),
        overwrite = TRUE
      )
      writeRaster(
        rast_undesirable_dev_rescaled,
        file.path(
          robustness_dir,
          paste0(
            "Undesirable_deviation_sum_of_change_",
            group_name,
            "_rescaled.tif"
          )
        ),
        overwrite = TRUE
      )
    }

    cat("Finished processing ES group:", group_name, "\n")
  }

  future::plan(future::sequential)
}

#' create_robustness_maps
#' Create bivariate robustness maps with different scaling and classification methods
#'
#' This function creates bivariate maps showing performance vs variation using different
#' scaling methods (min-max, quantile-based, robust MAD, winsorized) and classification
#' methods (equal intervals, quantile breaks, Fisher-Jenks). It saves maps, legends,
#' and palette information as JSON files.
#'
#' @param robustness_dir Character: Directory containing the robustness rasters and where outputs will be saved
#' @param group_names Character vector: Names of the groups to process (e.g., "ES", "BD", "ES_and_BD")
#' @param performance_file_pattern Character: Pattern for performance raster files (default: "Mean_sum_of_change_")
#' @param variation_file_pattern Character: Pattern for variation raster files (default: "Undesirable_deviation_sum_of_change_")
#' @param scaling_methods List: Named list of scaling methods with var_col, perf_col, and description
#' @param classification_methods List: Named list of classification methods with style and description
#' @param palette Character: Name of the biscale palette to use (default: "BlueGold")
#' @param num_classes Integer: Number of classes for bivariate classification (default: 3)
#' @param map_width Numeric: Width of output maps in cm (default: 25)
#' @param map_height Numeric: Height of output maps in cm (default: 20)
#' @param map_dpi Numeric: Resolution of output maps (default: 300)
#' @param save_json Logical: Whether to save palette information as JSON (default: TRUE)
#' @param quantile_probs Numeric vector: Quantiles to use for quantile-based scaling (default: c(0.02, 0.98))
#' @param robust_probs Numeric vector: Quantiles to use for robust scaling rescaling (default: c(0.01, 0.99))
#' @param winsor_probs Numeric vector: Quantiles to use for winsorization (default: c(0.05, 0.95))
#' @param ProjCH Character: Coordinate reference system to project rasters to (default: bash_vars$ProjCH)
#' @param verbose Logical: Whether to print progress messages (default: TRUE)
#'
#' @return NULL (saves files to disk)
#' @export
create_robustness_maps <- function(
  robustness_dir,
  group_names = c("ES"),
  performance_file_pattern = "Mean_sum_of_change_",
  variation_file_pattern = "Undesirable_deviation_sum_of_change_",
  scaling_methods = list(
    "minmax" = list(
      var_col = "Var_norm",
      perf_col = "Performance",
      description = "Min-Max Scaling"
    ),
    "quantile" = list(
      var_col = "Var_norm_quantile",
      perf_col = "Performance_norm_quantile",
      description = "Quantile-based Scaling"
    ),
    "robust" = list(
      var_col = "Var_norm_robust",
      perf_col = "Performance_norm_robust",
      description = "Robust (MAD) Scaling"
    ),
    "winsor" = list(
      var_col = "Var_norm_winsor",
      perf_col = "Performance_norm_winsor",
      description = "Winsorized Scaling"
    )
  ),
  classification_methods = list(
    "equal" = list(
      style = "equal",
      description = "Equal Intervals"
    ),
    "quantile" = list(
      style = "quantile",
      description = "Quantile Breaks"
    ),
    "fisher" = list(
      style = "fisher",
      description = "Fisher-Jenks"
    )
  ),
  palette = "BlueGold",
  num_classes = 3,
  map_width = 25,
  map_height = 20,
  map_dpi = 300,
  save_json = TRUE,
  quantile_probs = c(0.02, 0.98),
  robust_probs = c(0.01, 0.99),
  winsor_probs = c(0.05, 0.95),
  ProjCH = config$ProjCH,
  verbose = TRUE
) {
  # Ensure robustness directory exists
  ensure_dir(robustness_dir)

  for (group in group_names) {
    if (verbose) {
      cat(paste0("Processing group: ", group, "\n"))
    }

    # Load rasters
    performance_file <- file.path(
      robustness_dir,
      paste0(performance_file_pattern, group, ".tif")
    )
    variation_file <- file.path(
      robustness_dir,
      paste0(variation_file_pattern, group, ".tif")
    )

    # Check if files exist
    if (!file.exists(performance_file)) {
      warning(paste("Performance file not found:", performance_file))
      next
    }
    if (!file.exists(variation_file)) {
      warning(paste("Variation file not found:", variation_file))
      next
    }

    Performance <- rast(performance_file)
    Variation <- rast(variation_file)

    # project to desired CRS if needed
    if (!is.null(ProjCH)) {
      Performance <- project(Performance, ProjCH)
      Variation <- project(Variation, ProjCH)
    }

    # Stack the Performance and undesirable deviation rasters
    Combined_var <- c(Performance, Variation)
    names(Combined_var) <- c("Performance", "Var")

    # Get raster extent for consistent plotting
    raster_ext <- ext(Performance)
    x_range <- raster_ext$xmax - raster_ext$xmin
    y_range <- raster_ext$ymax - raster_ext$ymin
    data_aspect_ratio <- y_range / x_range

    # Convert to df
    Mean_var_df <- as.data.frame(Combined_var, xy = TRUE)

    # Apply all scaling methods

    # 1. Quantile-based rescaling
    var_p_low <- quantile(Mean_var_df$Var, quantile_probs[1], na.rm = TRUE)
    var_p_high <- quantile(Mean_var_df$Var, quantile_probs[2], na.rm = TRUE)
    Mean_var_df$Var_norm_quantile <- pmax(
      0,
      pmin(1, (Mean_var_df$Var - var_p_low) / (var_p_high - var_p_low))
    )

    # For Performance (-1 to 1 range, preserving zero)
    perf_pos <- Mean_var_df$Performance[Mean_var_df$Performance >= 0]
    perf_neg <- Mean_var_df$Performance[Mean_var_df$Performance < 0]

    if (length(perf_pos) > 0) {
      perf_pos_p_high <- quantile(perf_pos, quantile_probs[2], na.rm = TRUE)
    } else {
      perf_pos_p_high <- 0
    }

    if (length(perf_neg) > 0) {
      perf_neg_p_low <- quantile(perf_neg, quantile_probs[1], na.rm = TRUE)
    } else {
      perf_neg_p_low <- 0
    }

    Mean_var_df$Performance_norm_quantile <- ifelse(
      Mean_var_df$Performance >= 0,
      pmin(1, Mean_var_df$Performance / perf_pos_p_high),
      pmax(-1, Mean_var_df$Performance / abs(perf_neg_p_low))
    )

    # 2. Robust scaling using median and MAD
    var_median <- median(Mean_var_df$Var, na.rm = TRUE)
    var_mad <- mad(Mean_var_df$Var, na.rm = TRUE)
    var_robust_scaled <- (Mean_var_df$Var - var_median) / var_mad

    var_robust_min <- quantile(var_robust_scaled, robust_probs[1], na.rm = TRUE)
    var_robust_max <- quantile(var_robust_scaled, robust_probs[2], na.rm = TRUE)
    Mean_var_df$Var_norm_robust <- pmax(
      0,
      pmin(
        1,
        (var_robust_scaled - var_robust_min) / (var_robust_max - var_robust_min)
      )
    )

    perf_median <- median(Mean_var_df$Performance, na.rm = TRUE)
    perf_mad <- mad(Mean_var_df$Performance, na.rm = TRUE)
    perf_robust_scaled <- (Mean_var_df$Performance - perf_median) / perf_mad

    perf_robust_pos <- perf_robust_scaled[perf_robust_scaled >= 0]
    perf_robust_neg <- perf_robust_scaled[perf_robust_scaled < 0]

    if (length(perf_robust_pos) > 0) {
      perf_robust_pos_max <- quantile(
        perf_robust_pos,
        robust_probs[2],
        na.rm = TRUE
      )
    } else {
      perf_robust_pos_max <- 1
    }

    if (length(perf_robust_neg) > 0) {
      perf_robust_neg_min <- quantile(
        perf_robust_neg,
        robust_probs[1],
        na.rm = TRUE
      )
    } else {
      perf_robust_neg_min <- -1
    }

    Mean_var_df$Performance_norm_robust <- ifelse(
      perf_robust_scaled >= 0,
      pmin(1, perf_robust_scaled / perf_robust_pos_max),
      pmax(-1, perf_robust_scaled / abs(perf_robust_neg_min))
    )

    # 3. Winsorized scaling
    var_p_low_w <- quantile(Mean_var_df$Var, winsor_probs[1], na.rm = TRUE)
    var_p_high_w <- quantile(Mean_var_df$Var, winsor_probs[2], na.rm = TRUE)
    var_winsorized <- pmax(var_p_low_w, pmin(var_p_high_w, Mean_var_df$Var))

    Mean_var_df$Var_norm_winsor <- (var_winsorized -
      min(var_winsorized, na.rm = TRUE)) /
      (max(var_winsorized, na.rm = TRUE) - min(var_winsorized, na.rm = TRUE))

    perf_p_low_w <- quantile(
      Mean_var_df$Performance,
      winsor_probs[1],
      na.rm = TRUE
    )
    perf_p_high_w <- quantile(
      Mean_var_df$Performance,
      winsor_probs[2],
      na.rm = TRUE
    )
    perf_winsorized <- pmax(
      perf_p_low_w,
      pmin(perf_p_high_w, Mean_var_df$Performance)
    )

    Mean_var_df$Performance_norm_winsor <- ifelse(
      perf_winsorized >= 0,
      perf_winsorized / max(perf_winsorized, na.rm = TRUE),
      -1 * (perf_winsorized / min(perf_winsorized, na.rm = TRUE))
    )

    # 4. Min-max scaling (original)
    Mean_var_df$Var_norm <- (Mean_var_df$Var -
      min(Mean_var_df$Var, na.rm = TRUE)) /
      (max(Mean_var_df$Var, na.rm = TRUE) - min(Mean_var_df$Var, na.rm = TRUE))

    Mean_var_df$Performance <- ifelse(
      Mean_var_df$Performance >= 0,
      Mean_var_df$Performance / max(Mean_var_df$Performance, na.rm = TRUE),
      -1 *
        (Mean_var_df$Performance / min(Mean_var_df$Performance, na.rm = TRUE))
    )

    # Loop over scaling and classification methods
    for (method_name in names(scaling_methods)) {
      method <- scaling_methods[[method_name]]

      for (class_name in names(classification_methods)) {
        classification <- classification_methods[[class_name]]

        if (verbose) {
          cat(paste0(
            "Processing ",
            group,
            " with ",
            method$description,
            " and ",
            classification$description,
            "\n"
          ))
        }

        # Check if required columns exist
        if (!method$perf_col %in% names(Mean_var_df)) {
          warning(paste(
            "Performance column",
            method$perf_col,
            "not found. Skipping",
            method_name
          ))
          next
        }

        if (!method$var_col %in% names(Mean_var_df)) {
          warning(paste(
            "Variance column",
            method$var_col,
            "not found. Skipping",
            method_name
          ))
          next
        }

        # Create temporary dataframe with standard column names
        temp_df <- Mean_var_df
        temp_df$x_var <- Mean_var_df[[method$perf_col]]
        temp_df$y_var <- Mean_var_df[[method$var_col]]
        temp_df$x <- Mean_var_df$x
        temp_df$y <- Mean_var_df$y

        # Bivariate classification
        Mean_var_class <- bi_class(
          temp_df,
          x = x_var,
          y = y_var,
          style = classification$style,
          dim = num_classes
        )

        # Create breaks
        Mean_var_breaks <- bi_class_breaks(
          temp_df,
          x = x_var,
          y = y_var,
          style = classification$style,
          dim = num_classes,
          dig_lab = 2,
          split = FALSE
        )

        # Create map with ggplot2
        Full_robustness_map <- ggplot() +
          geom_raster(
            data = Mean_var_class,
            aes(x = x, y = y, fill = bi_class),
            show.legend = FALSE
          ) +
          bi_scale_fill(
            pal = palette,
            dim = num_classes,
            flip_axes = FALSE,
            rotate_pal = FALSE,
            na.value = "white"
          ) +
          coord_equal() +
          theme_void()

        # Save map using ggsave with same dimensions as base R
        map_filename <- file.path(
          robustness_dir,
          paste0(
            group,
            "_",
            method_name,
            "_",
            class_name,
            "_robustness_map.png"
          )
        )

        ggsave(
          plot = Full_robustness_map,
          filename = map_filename,
          width = map_width,
          height = map_height,
          units = "cm",
          dpi = map_dpi,
          bg = "transparent"
        )

        # Create and save legend
        Full_map_legend <- bi_legend(
          pal = palette,
          flip_axes = FALSE,
          rotate_pal = FALSE,
          dim = num_classes,
          breaks = Mean_var_breaks,
          xlab = "Mean sum of change",
          ylab = "Norm. undesirable deviation",
          arrow = FALSE
        ) +
          theme(
            text = element_text(size = 8, family = "FiraSans"),
            axis.text.x = element_text(angle = -25, hjust = 0),
            axis.text.y = element_text(angle = -25, hjust = 0.25),
            panel.background = element_rect(fill = 'transparent'),
            plot.background = element_rect(fill = 'transparent', color = NA)
          )

        ggsave(
          Full_map_legend,
          filename = file.path(
            robustness_dir,
            paste0(
              group,
              "_",
              method_name,
              "_",
              class_name,
              "_robustness_legend.png"
            )
          ),
          width = map_width,
          height = map_height,
          units = "cm",
          dpi = map_dpi,
          bg = "transparent"
        )

        # Get bivariate color palette for JSON
        bi_colors_plot <- bi_pal(
          pal = palette,
          dim = num_classes,
          flip_axes = FALSE,
          rotate_pal = FALSE
        )

        bi_colors <- ggplot_build(bi_colors_plot)$data[[1]]$fill

        bi_class_values <- paste0(
          rep(1:num_classes, each = num_classes),
          "-",
          rep(1:num_classes, num_classes)
        )

        if (length(bi_colors) != length(bi_class_values)) {
          bi_colors <- biscale:::bi_pal_manual(
            pal = palette,
            dim = num_classes,
            flip_axes = FALSE,
            rotate_pal = FALSE
          )
        }

        # Save JSON palette information
        if (save_json) {
          palette_info <- list()

          for (i in 1:length(bi_class_values)) {
            class_id <- bi_class_values[i]
            color_code <- bi_colors[i]

            parts <- strsplit(class_id, "-")[[1]]
            x_pos <- as.numeric(parts[1])
            y_pos <- as.numeric(parts[2])

            perf_labels <- c("Low", "Medium", "High")
            var_labels <- c("Low", "Medium", "High")

            interpretation <- paste0(
              perf_labels[x_pos],
              " Performance, ",
              var_labels[y_pos],
              " Variation"
            )

            palette_info[[class_id]] <- list(
              color = color_code,
              interpretation = interpretation
            )
          }

          palette_json <- toJSON(palette_info, pretty = TRUE, auto_unbox = TRUE)

          json_filename <- file.path(
            robustness_dir,
            paste0(
              group,
              "_",
              method_name,
              "_",
              class_name,
              "_palette_info.json"
            )
          )

          write(palette_json, file = json_filename)

          if (verbose) {
            cat("Saved palette information to:", json_filename, "\n")
          }
        }
      }
    }

    if (verbose) {
      cat(paste0("Completed processing for group: ", group, "\n"))
    }
  }

  invisible(NULL)
}


#' summarize_ES_outputs:
#' This function normalises the ES layers using the global minmax values,
#' calculates summary statistics for the rescaled layers,
#' and calculates the areas of land that fall within discrete classes of the ES provision values.
#' It saves the results to the specified directories.
#' @param metrics Vector of metrics to calculate must use built-in function arguments from terra::global
#' i.e. "max", "min", "mean", "sum", "range", "rms" (root mean square), "sd",
#'  "std" (population sd, using n rather than n-1),
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param ESs_to_summarise Vector of ESs to summarise
#' @param ES_rescaling_dir Directory to load the global minmaxs values used for normalization
#' @param ES_summary_stats_dir Directory to save the summary statistics for individual ESs
#' @param ES_summarisation_dir Directory to save the summary statistics for all ESs
#' @param recalc_minmax = FALSE,
#' @param recalc_rescaled_layers = FALSE,
#' @param recalc_summary_stats = TRUE,
#' @param recalc_perc_area = FALSE,
#' @param recalc_es_chg_maps = TRUE,
#' @param recalc_norm_chg_maps = TRUE,
#' @param recalc_chg_perc_area = TRUE,
#' @param recalc_cumulative_summary = FALSE, # NEW PARAMETER
#' @param Sim_ctrl_tbl_path Character: The path to the simulation control table
#' @param Summarize_across_ES Logical: If TRUE, the function will summarize the results across all ESs
#' @param Rescale_cross_ES_results Logical: If TRUE, the function will rescale the summary results across all ESs
#' @param parallel Logical: If TRUE, the function will run in parallel using the number of workers specified in the config file
#' @param mask Character: The path to the mask to apply to the layers
#' @param sample_size Integer: The number of samples to take from each layer for calculating global breaks
#' @param perc_area_break_types Character: The types of breaks to calculate for classified area, options are "quantile", "fisher" and "regular"
#' @param save_provision_area_plots Logical: Whether to save the classified area plots
#' @param save_break_samples Logical: Whether to save the samples used for calculating breaks
#' @return: NULL
#'
#' @export
#
summarize_ES_outputs <- function(
  ProjCH = "EPSG:2056",
  metrics = c("sum", "mean", "sd"),
  ES_layer_paths,
  ESs_to_summarise,
  ES_rescaling_dir,
  ES_summary_stats_dir,
  ES_summarisation_dir,
  recalc_minmax = FALSE,
  recalc_rescaled_layers = FALSE,
  recalc_summary_stats = TRUE,
  recalc_perc_area = FALSE,
  recalc_es_chg_maps = TRUE,
  recalc_norm_chg_maps = TRUE,
  recalc_chg_perc_area = TRUE,
  recalc_cumulative_summary = FALSE, # NEW PARAMETER
  Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
  Summarize_across_ES = TRUE,
  Rescale_cross_ES_results = TRUE,
  parallel = TRUE,
  n_workers = NULL,
  mask,
  sample_size = 50000,
  perc_area_break_types = c("quantile", "fisher", "regular"),
  perc_area_focus_regions = list(
    "All" = NULL,
    "High_prov" = c("High_prov"),
    "Low_prov" = c("Low_prov")
  ),
  save_provision_area_plots = TRUE,
  save_break_samples = FALSE
) {
  # Set parallel processing options
  if (parallel == FALSE) {
    future::plan(future::sequential)
  } else {
    # Use provided n_workers or default to all available cores
    if (is.null(n_workers) || n_workers == 0) {
      n_workers <- future::availableCores()
    }
    future::plan(future::multisession, workers = n_workers)
  }

  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)

  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "EI",
    "EI_"
  )

  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "GREX",
    "GR_EX"
  )

  Mask_summary_stats_dir <- file.path(ES_summary_stats_dir, names(mask))

  # if the directory does not exist, create it
  ensure_dir(Mask_summary_stats_dir)

  # Loop over the ESs specified and calculate the summary metrics specified for the
  # rescaled values for each time step under all configurations
  # for each configuration and time-step
  lapply(ESs_to_summarise, function(ES) {
    # convert ES to lowercase
    ES <- tolower(ES)

    # subset to ES records
    ES_layer_paths <- ES_layer_paths[ES_layer_paths$ES == ES, ]

    # Step 1: Calculate minmaxs for this ES
    # Create the directory to save the minmax values for each ES
    Mask_rescaling_dir <- file.path(ES_rescaling_dir, names(mask))
    ensure_dir(Mask_rescaling_dir)

    # create an output file path for this ES
    ES_minmax_path <- file.path(Mask_rescaling_dir, paste0(ES, "_minmaxs.rds"))
    cat("Step 1: Calculating minmaxs for all configs under", ES, "\n")
    if (file.exists(ES_minmax_path) & recalc_minmax == FALSE) {
      cat(
        "a previous minmax file exists and recalc_minmax == FALSE: Loading existing file ",
        ES_minmax_path,
        "\n"
      )
      ES_minmax <- readRDS(ES_minmax_path)
    } else if (!file.exists(ES_minmax_path) | recalc_minmax == TRUE) {
      cat(
        "recalc_minmax == TRUE or no previous file exists: Calculating minmaxs"
      )
      ES_minmax <- calc_minmaxs(
        ES_layer_paths = ES_layer_paths,
        save_path = ES_minmax_path,
        report_NAs = FALSE,
        mask = mask,
        map_path = "Path",
        ES = ES
      )
    }

    # Step 2: Calculate the global min max values for this ES
    cat("Step 2: Summarising global minmaxs for", ES, "\n")
    ES_global_minmaxs <- calc_global_minmaxs(ES_minmax = ES_minmax)

    # Step 3. Rescale the ES layers
    cat("Step 3: Performing min-max rescaling of layers for", ES, "\n")
    if (recalc_rescaled_layers == TRUE) {
      cat("recalc_rescaled_layers == TRUE: Recaling all layers")
      # Apply the normalisation function to the ES layers
      normalise_layers_PATCH(
        ES_layer_paths = ES_layer_paths,
        mask = mask,
        ES_global_minmaxs = ES_global_minmaxs,
        ProjCH = ProjCH,
        input_path_name = "Path",
        raster_path_name = "tif_path",
        image_path_name = "png_path",
        overwrite = TRUE,
        low_color = "#007CDC",
        high_color = "#FFEA2A",
        mid_color = NULL,
        low_value = 0,
        high_value = 1,
        mid_value = NULL
      )
    } else if (recalc_rescaled_layers == FALSE) {
      cat(
        "Recalculation of rescaled layers is set to FALSE, identifying which layers already exist for",
        ES,
        "\n"
      )

      # Count how many of the rescaled layers for this ES exist and how many
      # are missing
      ES_layer_paths$Exists <- sapply(ES_layer_paths$tif_path, file.exists)

      cat(
        "Number of rescaled raster layers for",
        ES,
        "that exist:",
        sum(ES_layer_paths$Exists),
        "\n"
      )
      cat(
        "Number of rescaled raster layers for",
        ES,
        "that are missing:",
        sum(!ES_layer_paths$Exists),
        "\n"
      )

      # Remove the existing layers from the dataframe
      ES_layer_paths_remaining <- ES_layer_paths[!ES_layer_paths$Exists, ]

      cat(
        "Existing layers have been removed from the list to be rescaled, performing rescaling of remaining layers \n"
      )

      # If there are no layers to rescale, skip to the next ES
      if (nrow(ES_layer_paths_remaining) == 0) {
        cat(
          "Recalculation of rescaled layers is set to FALSE and all layers already exist for",
          ES,
          "skipping to next step \n"
        )
      } else if (nrow(ES_layer_paths_remaining) > 0) {
        cat(
          "Recalculation of rescaled layers is set to FALSE, rescaling remaining layers for",
          ES,
          "\n"
        )

        # If there are layers to rescale, apply the normalisation function
        normalise_layers_PATCH(
          ES_layer_paths = ES_layer_paths_remaining,
          mask = mask,
          ES_global_minmaxs = ES_global_minmaxs,
          ProjCH = ProjCH,
          input_path_name = "Path",
          raster_path_name = "tif_path",
          image_path_name = "png_path",
          overwrite = TRUE,
          low_color = "#007CDC",
          high_color = "#FFEA2A",
          mid_color = NULL,
          low_value = 0,
          high_value = 1,
          mid_value = NULL
        )
      }
    }

    # Step 4: Calculate summary statistics for the rescaled layers
    # create path to save the combined results for this ES
    cat("Step 4: Calculating summary stats for", ES, "\n")
    ES_sum_stats_path <- file.path(
      Mask_summary_stats_dir,
      paste0(tolower(ES), "_summary_stats.rds")
    )

    # check if the summary file already exists & if recalc_summary_stats == FALSE then skip
    # otherwise if the file does not exist or recalc_summary_stats = TRUE then re-calculate
    if (file.exists(ES_sum_stats_path) & recalc_summary_stats == FALSE) {
      cat("Summary stats file already exists for", ES, "\n")
      cat("Recalculation of summary is set to FALSE, skipping to next step \n")
    } else if (!file.exists(ES_sum_stats_path) | recalc_summary_stats == TRUE) {
      if (!file.exists(ES_sum_stats_path)) {
        cat(
          "Summary stats file does not exist for",
          ES,
          ": performing recalculation",
          "\n"
        )
      }
      if (recalc_summary_stats == TRUE) {
        cat(
          "Recalculation of summary is set to TRUE, calculating Summary stats for",
          ES,
          "\n"
        )
      }

      # Apply function to calculate summary statistics for the rescaled layers
      ES_sum_stats <- calc_summary_stats(ES_layer_paths, metrics = metrics)
      # Save the results for the current ES
      saveRDS(ES_sum_stats, file = ES_sum_stats_path)
      cat(
        "Summary stats for",
        ES,
        "have been saved to",
        ES_sum_stats_path,
        "\n"
      )
    }

    # Step 5: Classify the rescaled layers and calculate the area of each class
    cat(
      "Step 5: Calculating % area of land in discrete classes of ES provision for",
      ES,
      "\n"
    )
    if (recalc_perc_area == TRUE) {
      # Calculate the areas in classes for the rescaled layers
      Calculate_areas_in_classes(
        ES_layer_paths = ES_layer_paths,
        sample_size = sample_size,
        save_samples = save_break_samples,
        break_types = perc_area_break_types,
        focus_regions = perc_area_focus_regions,
        save_plots = save_provision_area_plots,
        parallel = parallel,
        n_workers = n_workers,
        input_map_path_name = "tif_path",
        output_tbl_path_name = "perc_area_path",
        regular_min = 0,
        regular_max = 1,
        regular_by = 0.1
      )
      cat(
        "% areas in classes for",
        ES,
        "have been calculated and saved to JSON files \n"
      )
    }

    # Step 6: Produce rasters of change in ES provision over time
    # and calculate minmax values
    cat(
      "Step 6: Creating maps of change in provision and calculating min max values for",
      ES,
      "\n"
    )

    ES_chg_minmax_path <- file.path(
      Mask_rescaling_dir,
      paste0(ES, "_chg_minmaxs.rds")
    )

    if (file.exists(ES_chg_minmax_path) & recalc_es_chg_maps == FALSE) {
      cat(
        "recalc_es_chg_maps == FALSE and a a file containing the min max values of change in provision exists: Loading existing file \n"
      )
      ES_chg_minmaxs <- readRDS(ES_chg_minmax_path)
    } else if (!file.exists(ES_chg_minmax_path) | recalc_es_chg_maps == TRUE) {
      cat(
        "Either a minmax file for change in provision does not exist or recalc_es_chg_maps == TRUE: Producing maps of change in ES provision over time for",
        ES,
        "\n"
      )
      ES_chg_minmaxs <- Produce_change_map(
        ES_layer_paths = ES_layer_paths,
        diff_vs_init_path_name = "chg_tif_init_path",
        sum_chg_path_name = "sum_chg_tif_path",
        ES_chg_minmax_path = ES_chg_minmax_path
      )
    }

    # seperate the sum_change minmaxs from the time step minmaxs
    ES_chg_minmaxs <- ES_chg_minmaxs[names(ES_chg_minmaxs) != "sum_change"]

    # Step 7: Summarize global min max values of change in provision vs initial time step
    cat(
      "Step 7: Summarising global min maxs of change in provision for",
      ES,
      "\n"
    )
    ES_chg_global_minmaxs <- calc_global_minmaxs(ES_minmax = ES_chg_minmaxs)

    # remove all entries with Time_step == 2020 because change layers do not exist for these
    Chg_layer_paths <- ES_layer_paths[ES_layer_paths$Time_step != 2020, ]

    # Step 8: rescale rasters of change in ES provision
    cat("Step 8: Rescaling maps of ES change in provision for", ES, "")
    if (recalc_norm_chg_maps == TRUE) {
      normalise_layers(
        ES_layer_paths = Chg_layer_paths,
        mask = mask,
        ES_global_minmaxs = ES_chg_global_minmaxs,
        ProjCH = ProjCH,
        input_path_name = "chg_tif_init_path",
        raster_path_name = "chg_tif_path",
        image_path_name = "chg_png_path",
        overwrite = TRUE,
        low_color = "#D10062",
        high_color = "#D1FF61",
        mid_color = "#B2B2B2",
        low_value = -1,
        high_value = 1,
        mid_value = 0
      )
    }

    # Step 9: Classify the change layers and calculate the area of each class
    cat(
      "Step 9: Calculating % area of land in discrete classes of change in ES provision for",
      ES,
      "\n"
    )
    if (recalc_chg_perc_area == TRUE) {
      # Calculate the areas in classes for the rescaled layers
      Calculate_areas_in_classes(
        ES_layer_paths = Chg_layer_paths,
        sample_size = sample_size,
        save_samples = save_break_samples,
        break_types = perc_area_break_types,
        focus_regions = perc_area_focus_regions,
        save_plots = save_provision_area_plots,
        parallel = parallel,
        n_workers = n_workers,
        input_map_path_name = "chg_tif_path",
        output_tbl_path_name = "area_chg_path",
        regular_min = -1,
        regular_max = 1,
        regular_by = 0.25
      )
      cat(
        "Areas in classes for",
        ES,
        "have been calculated and saved to JSON files \n"
      )
    }
  }) # close loop over ESs_to_summarise

  # Summarize cumulative ES change maps
  if (recalc_cumulative_summary == TRUE) {
    cat("Step 10: Summarizing cumulative ES change maps\n")

    # Define ES groups - you can customize this or make it a parameter
    ES_groups <- list("ES" = ESs_to_summarise)

    summarize_cumulative_ES_change_maps(
      sum_chg_rast_dir = file.path(config$OutputDir, config$SumChgRasterDir),
      robustness_dir = file.path(config$OutputDir, config$RobustnessDir),
      ES_groups = ES_groups,
      spatial_summary_metrics = c(
        #"mean",
        "undesirable_deviation"
      )
    )

    cat(
      "Cumulative ES change maps summarized and saved to",
      file.path(config$OutputDir, config$RobustnessDir),
      "\n"
    )
  }

  # if parallel == TRUE set back to sequential processing
  if (parallel == TRUE) {
    future::plan(future::sequential)
  }

  # Summarise the results across all ESs
  if (Summarize_across_ES == TRUE) {
    # Produce a single df summarising the results across all ESs
    ES_sum_stats <- lapply(ESs_to_summarise, function(ES) {
      # Load the summary stats for the current ES
      ES_sum_stats <- readRDS(file.path(
        Mask_summary_stats_dir,
        paste0(tolower(ES), "_summary_stats.rds")
      ))

      # Bind the list of dataframes for each time point together
      ES_sum_stats <- do.call(rbind, ES_sum_stats)

      # remove the exists column
      ES_sum_stats$Exists <- NULL

      return(ES_sum_stats)
    })

    #rbind list of dataframes into a single dataframe
    ES_sum_stats <- do.call(rbind, ES_sum_stats)

    # Use the simulation control table to add the scenario_ID to the ES_sum_stats
    ES_sum_stats$Scenario <- sapply(ES_sum_stats$Config_ID, function(x) {
      Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == x, "Scenario_ID.string"]
    })

    # Convert metric columns to numeric
    ES_sum_stats[, metrics] <- lapply(ES_sum_stats[, metrics], as.numeric)

    # create dir for saving summarisation results
    Mask_summarisation_dir <- file.path(ES_summarisation_dir, names(mask))
    if (!dir.exists(Mask_summarisation_dir)) {
      dir.create(Mask_summarisation_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Save the ES_sum_stats
    saveRDS(
      ES_sum_stats,
      file = file.path(Mask_summarisation_dir, "ES_summary_stats.rds")
    )

    cat("Summary stats across all ES have been calculated and saved \n")

    # if Rescale_results is TRUE, then normalise the results and save a copy
    if (Rescale_cross_ES_results == TRUE) {
      cat(
        "rescaled_summary_stats is set to TRUE, normalising the summary stats \n"
      )

      # Normalise the summary stats
      ES_sum_stats_norm <- lapply(ESs_to_summarise, function(ES) {
        # Get summary values for current ES
        ES_sum_stats_ES <- ES_sum_stats[ES_sum_stats$ES == ES, ]

        # Loop over the metrics
        for (metric in metrics) {
          # Initialize min and max values
          ES_min <- min(ES_sum_stats_ES[[metric]], na.rm = TRUE)
          ES_max <- max(ES_sum_stats_ES[[metric]], na.rm = TRUE)

          # Rescale the summary values
          ES_sum_stats_ES[[paste0(metric, "_rescale")]] <- (ES_sum_stats_ES[[
            metric
          ]] -
            ES_min) /
            (ES_max - ES_min)

          # Remove the original summary values
          ES_sum_stats_ES[[metric]] <- NULL

          # Rename the rescaled summary values
          names(ES_sum_stats_ES)[
            names(ES_sum_stats_ES) == paste0(metric, "_rescale")
          ] <- metric
        }

        return(ES_sum_stats_ES)
      })
      names(ES_sum_stats_norm) <- ESs_to_summarise

      # bind the list of dataframes into a single dataframe
      ES_sum_stats_norm <- as.data.frame(rbindlist(ES_sum_stats_norm))

      # Save the rescaled summary stats
      saveRDS(
        ES_sum_stats_norm,
        file = file.path(
          Mask_summarisation_dir,
          "ES_summary_stats_rescaled.rds"
        )
      )
    }
  }

  cat("Finished normalising and summarising the ES layers \n")
}


#' Summarize_for_masks
#' Upper level function to summarise the ESs for a list of supplied masks
#' #' @param config List: Configuration list containing the ESs to summarise,
#' #' ES nesting, ES file names, and other parameters
#' #' @param output_dir Character: The directory where the web platform files are stored
#' #' @param ProjCH string of CRS
#' #' @param ES_rescaling_dir Character: The directory where the rescaled ES layers are stored
#' #' @param masks List: A named list of masks to apply to the ES layers if the
#' #  name is full_extent then the list item can be empty and the summary will be performed on the entire map
#' #' @param parallel Logical: If TRUE, parallel processing will be used (read from config)
#' #' @param n_workers Integer: Number of workers for parallel processing (read from config)
Summarise_for_masks <- function(
  config = config,
  ProjCH,
  ES_rescaling_dir = ES_rescaling_dir,
  ES_summary_stats_dir = ES_summary_stats_dir,
  ES_summarisation_dir = ES_summarisation_dir,
  masks = list("canton" = file.path(Mask_dir, "Canton_mask.shp")),
  # Add arguments that will be passed to summarize_ES_outputs
  metrics = c("sum", "mean", "sd"),
  recalc_minmax = FALSE,
  recalc_rescaled_layers = FALSE,
  recalc_summary_stats = FALSE,
  recalc_perc_area = TRUE,
  recalc_es_chg_maps = FALSE,
  recalc_norm_chg_maps = FALSE,
  recalc_chg_perc_area = TRUE,
  recalc_cumulative_summary = FALSE,
  Summarize_across_ES = FALSE,
  Rescale_cross_ES_results = FALSE,
  parallel = NULL,
  n_workers = NULL
) {
  # Set parallel processing from config if not explicitly provided
  if (is.null(parallel)) {
    parallel <- ifelse(
      is.null(config$Parallel) || config$Parallel == "",
      FALSE,
      config$Parallel
    )
  }

  # Set number of workers from config if not explicitly provided
  if (is.null(n_workers)) {
    n_workers <- ifelse(
      is.null(config$NWorkers) || config$NWorkers == "" || config$NWorkers == 0,
      future::availableCores(),
      config$NWorkers
    )
  }

  cat(paste0(
    "Parallel processing: ",
    parallel,
    " (workers: ",
    n_workers,
    ")\n"
  ))
  # loop over the masks
  for (i in seq_along(masks)) {
    cat(paste0("Processing ES outputs for mask: ", names(masks)[i], "\n"))

    # apply function to prepare DF
    ES_layer_paths <- prepare_analysis_df(
      input_dir = config$InputDir,
      image_dir = config$ImageDir,
      chg_image_dir = config$ChgImageDir,
      raster_dir = config$RasterDir,
      chg_raster_dir = config$ChgRasterDir,
      chg_raster_init = config$ChgRasterInit,
      sum_chg_raster_dir = config$SumChgRasterDir,
      perc_area_dir = config$PercAreaDir,
      area_chg_dir = config$AreaChgDir,
      sample_dir = "sample_data",
      base_dir = config$OutputDir,
      Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
      ProjCH = ProjCH,
      ESs_to_summarise = config$ESs_to_summarise,
      ES_nesting = config$ES_nesting,
      ES_file_names = config$ES_file_names,
      mask = masks[i]
    )

    # Add output dir to focus areas
    perc_area_focus_regions <- lapply(config$PercAreaFocusRegions, function(x) file.path(config$OutputDir, x))
    names(perc_area_focus_regions) <- names(config$PercAreaFocusRegions)

    # Normalise and summarise the ESs for this mask
    summarize_ES_outputs(
      ProjCH = ProjCH,
      metrics = metrics,
      ES_layer_paths = ES_layer_paths,
      ESs_to_summarise = config$ESs_to_summarise,
      ES_rescaling_dir = ES_rescaling_dir,
      ES_summary_stats_dir = ES_summary_stats_dir,
      ES_summarisation_dir = ES_summarisation_dir,
      recalc_minmax = recalc_minmax,
      recalc_rescaled_layers = recalc_rescaled_layers,
      recalc_summary_stats = recalc_summary_stats,
      recalc_perc_area = recalc_perc_area,
      recalc_es_chg_maps = recalc_es_chg_maps,
      recalc_norm_chg_maps = recalc_norm_chg_maps,
      recalc_chg_perc_area = recalc_chg_perc_area,
      recalc_cumulative_summary = recalc_cumulative_summary,
      Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
      Summarize_across_ES = Summarize_across_ES,
      Rescale_cross_ES_results = Rescale_cross_ES_results,
      parallel = parallel,
      n_workers = n_workers,
      sample_size = 50000,
      perc_area_break_types = c("regular"),
      perc_area_focus_regions = perc_area_focus_regions,
      save_provision_area_plots = TRUE,
      save_break_samples = TRUE,
      mask = masks[i]
    )
  }
}


### =========================================================================
### Main
### =========================================================================

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
config <- config$Summarisation # only summarisation variables

# if InputDir is not set, use bash variable
# $FUTURE_EI_OUTPUT_DIR/$ES_OUTPUT_BASE_DIR
if (is.null(config$InputDir) || config$InputDir == "") {
  config$InputDir <- file.path(
    bash_vars$FUTURE_EI_OUTPUT_DIR,
    bash_vars$ES_OUTPUT_BASE_DIR
  )
}
# if OutputDir is not set, use bash variable SUMMARISATION_OUTPUT_DIR
if (is.null(config$OutputDir) || config$OutputDir == "") {
  config$OutputDir <- file.path(
    bash_vars$FUTURE_EI_OUTPUT_DIR,
    bash_vars$SUMMARISATION_OUTPUT_DIR
  )
}
# Check if InputDir is now set
if (is.null(config$InputDir) || config$InputDir == "") {
  stop("InputDir nor ES_OUTPUT_BASE_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if InputDir exists
if (!dir.exists(config$InputDir)) {
  stop(cat("InputDir does not exist: ", config$InputDir, "\n"))
}
# Check if OutputDir is set
if (is.null(config$OutputDir) || config$OutputDir == "") {
  stop("OutputDir nor SUMMARISATION_OUTPUT_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if OutputDir exists and create if not
ensure_dir(config$OutputDir)

cat("Performing summarisation of ES outputs \n")

cat("Working directory is:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")

### Start by creating sub directories

# Mask directory
Mask_dir <- file.path(config$OutputDir, "Masks")
# if dir doesn't exist, create it
ensure_dir(Mask_dir)

# Dir for the processed ES layers
# Sub-dir for the summarisation results
ES_summarisation_dir <- file.path(config$OutputDir, "ES_summarisation")
ensure_dir(ES_summarisation_dir)

# Sub-dir for the normalisation results
ES_rescaling_dir <- file.path(ES_summarisation_dir, "ES_rescaling")
ensure_dir(ES_rescaling_dir)

# Sub-dir for the summary_stats results
ES_summary_stats_dir <- file.path(ES_summarisation_dir, "ES_summary_stats")
ensure_dir(ES_summary_stats_dir)

### =========================================================================
### Produce ES outputs for web platform
### =========================================================================

# run function to summarise for masks
Summarise_for_masks(
  config = config,
  ProjCH = config$ProjCH,
  ES_rescaling_dir = ES_rescaling_dir,
  ES_summary_stats_dir = ES_summary_stats_dir,
  ES_summarisation_dir = ES_summarisation_dir,
  masks = c(
    "canton" = file.path(Mask_dir, "Bern_Canton_w_lakes.shp")
  ),
  metrics = c("sum", "mean", "sd"),
  recalc_minmax = FALSE,
  recalc_rescaled_layers = FALSE,
  recalc_summary_stats = FALSE,
  recalc_perc_area = TRUE,
  recalc_es_chg_maps = TRUE,
  recalc_norm_chg_maps = TRUE,
  recalc_chg_perc_area = TRUE,
  recalc_cumulative_summary = TRUE,
  Summarize_across_ES = TRUE,
  Rescale_cross_ES_results = TRUE
)

### =========================================================================
### Robustness map
### =========================================================================

# # Define directory to save robustness maps
robustness_dir <- file.path(config$OutputDir, config$RobustnessDir)

# because the SDM outputs were calculated seperately and cumulative sums have not been calculated yet we need to do so here
sdm_chg_rast_dir <- "X:/CH_Kanton_Bern/03_Workspaces/02_SDM/sdm_summarisation/chg_raster_data"

# list all tif files that contain the string uzl_npa-all
sdm_files <- list.files(
  sdm_chg_rast_dir,
  pattern = "uzl_npa.*\\.tif$",
  full.names = TRUE
)

# remove the string sdm-XXXX- (where XXXX is a four digit numeric) from the start of the file names
unique_configs <- unique(gsub(
  "sdm-\\d{4}-",
  "",
  basename(sdm_files)
))

# path to the sum_chg_rast_dir
sum_chg_rast_dir <- file.path(
  config$OutputDir,
  config$SumChgRasterDir
)

# the extent of the sdm rasters does not match all others in one dimension, use one of the lulc rasters to modify the extent
es_test <- rast(
  "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/raster_data/hab-2030-rcp85-ref_peri_urban-ref-ei_soc-canton.tif"
)

# loop over the unique configs and calculate the cumulative sum of change rasters
for (config_name in unique_configs) {
  # list all files for this config
  config_files <- sdm_files[grepl(config_name, sdm_files)]

  # read in the rasters and calculate the cumulative sum
  rasters <- rast(config_files)
  cum_sum_raster <- sum(rasters)

  # resample to match the lulc_test raster
  cum_sum_raster <- resample(cum_sum_raster, es_test, method = "bilinear")

  # reproject to the project CRS
  cum_sum_raster <- project(cum_sum_raster, crs(es_test))

  # run camparegeom for the first ratser only
  if (config_name == unique_configs[1]) {
    compareGeom(cum_sum_raster, es_test, stopOnError = TRUE)
  }

  # define output file path
  output_file <- file.path(
    sum_chg_rast_dir,
    paste0("sdm-", config_name)
  )

  # write the cumulative sum raster to file
  terra::writeRaster(
    cum_sum_raster,
    filename = output_file,
    overwrite = TRUE
  )
}

# Apply function to summarize cumulative ES change maps
summarize_cumulative_ES_change_maps(
  sum_chg_rast_dir = file.path(config$OutputDir, config$SumChgRasterDir),
  robustness_dir = file.path(config$OutputDir, config$RobustnessDir),
  ES_groups = list(
    "BD-ES" = c('sdm', 'REC', 'CAR', 'NDR', 'POL', 'SDR', 'WY', 'FF', 'HAB')
  ),
  spatial_summary_metrics = c(
    "mean",
    #"stdev",
    #"mean-variance",
    "undesirable-deviation"
  )
)


# NOTE: FOR THE ES ONLY MAP I HAVE USED WINSOR NORMALISATION AND 'EQUAL' FOR THE MAP CATEGORIES

# # Create robustness maps with default settings
create_robustness_maps(
  robustness_dir = robustness_dir,
  group_names = c("BD", "BD-ES"), # , "BD", "BD-ES"
  verbose = TRUE,
  ProjCH = config$ProjCH,
  scaling_methods = list(
    "winsor" = list(
      var_col = "Var_norm_winsor",
      perf_col = "Performance_norm_winsor",
      description = "Winsorized Scaling"
    )
  ),
  classification_methods = list(
    "quantile" = list(
      style = "quantile",
      description = "Quantile Breaks"
    )
  )
)
