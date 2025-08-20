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
packs <- c("stringr", "terra", "future", "future.apply", "readxl",
           "data.table", "tidyr", "yaml", "dplyr", "viridis", "ggplot2",
           "tidyterra", "jsonlite", "magick", "grDevices", "classInt")
invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

### =========================================================================
### Functions
### =========================================================================

### Prepare a dataframe of information on the ES layers that are to be processed
#' Prepare a dataframe of information on the ES layers that are to be processed
#' Including the input data path for each config and time step, 
#' a path to save the normalized tif file, a path to save a png file of the rescaled layer
#' a path to save a JSON file of the areas of land that fall within discrete classes of the ES provision values
#' 
#' @param input_dir The directory where the input ES layers are stored
#' @param image_dir The directory where the output images will be saved
#' @param raster_dir The directory where the output rasters will be saved
#' @param classified_area_dir The directory where the classified area JSON files will be saved
#' @param base_dir The base directory where the output directories are located
#' @param Sim_ctrl_tbl_path The path to the simulation control table CSV file
#' @param ProjCH The projection CRS to use for the ES layers
#' @param ESs_to_summarise A character vector of ESs to summarise
#' @param ES_nesting A named list indicating whether each ES is nested (TRUE) or not (FALSE)
#' @param ES_file_names A named list of file names for each ES
#' @param map_masks A list of shapefiles to use as masks for the ES layers
#' @return A dataframe with information on the ES layers to be processed, including paths for input, output rasters, images, and classified areas
prepare_analysis_df <- function(
    input_dir = "F:/KB-outputs/ncp_output",
    image_dir = "map_images",
    raster_dir = "raster_data",
    classified_area_dir = "chart_data/classified_area",
    sample_dir = "sample_data",
    base_dir = web_platform_dir,
    Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
    ProjCH = ProjCH,
    ESs_to_summarise = config$ESs_to_summarise,
    ES_nesting = config$ES_nesting,
    ES_file_names = config$ES_file_names,
    mask
){
  
  # if the directories do not exist, create them
  if(!dir.exists(file.path(base_dir, image_dir))){
    dir.create(file.path(base_dir, image_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, raster_dir))){
    dir.create(file.path(base_dir, raster_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, classified_area_dir))){
    dir.create(file.path(base_dir, classified_area_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, sample_dir))){
    dir.create(file.path(base_dir, sample_dir), recursive = TRUE, showWarnings = FALSE)
  }

  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)
  
  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI_")
  
  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "GREX", "GR_EX")
  
    # convert all values in columns: Climate_scenario.string, Econ_scenario.string,
  # Pop_scenario.string and Scenario_ID.string to lower case
  Sim_ctrl_tbl$Climate_scenario.string <- tolower(Sim_ctrl_tbl$Climate_scenario.string)
  Sim_ctrl_tbl$Econ_scenario.string <- tolower(Sim_ctrl_tbl$Econ_scenario.string)
  Sim_ctrl_tbl$Pop_scenario.string <- tolower(Sim_ctrl_tbl$Pop_scenario.string)
  Sim_ctrl_tbl$Scenario_ID.string <- tolower(Sim_ctrl_tbl$Scenario_ID.string)
  
  # Get earliest scenario start date and latest end date
  Start_date <- min(Sim_ctrl_tbl$Scenario_start.real)
  End_date <- max(Sim_ctrl_tbl$Scenario_end.real)
  
  # Create seq of scenario time steps with Step_length.real
  Sim_time_steps <- seq(Start_date, End_date, by = Sim_ctrl_tbl$Step_length.real[1])
  
  # Vector IDs of configurations to be analysed
  # Note Manually adjust this to analyse specific configurations
  Config_IDs <- unique(Sim_ctrl_tbl$Simulation_num.)
  
  # Loop over ESs_to_summarise and create vectors of file paths for their layers
  # according to Config_IDs, Sim_time_steps and nesting
  ES_paths <- lapply(ESs_to_summarise, function(ES){
    
    # convert ES to lowercase
    ES_lower <- tolower(ES)
    
    # Loop over config_IDs returning paths as vector
    Config_paths <- lapply(Config_IDs, function(Config_ID){
      
      # use Config_ID to subset Sim_ctrl_tbl to get scenario details
      Config_details <- Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == Config_ID, ]
      
      # Create a vector of file paths for each time step
      Config_time_paths <- as.data.frame(rbindlist(lapply(Sim_time_steps, function(Time_step){
      
        # If ES is not nested then append time step to file path
        if(ES_nesting[[ES]] == FALSE){
          input_path <- file.path(input_dir, Config_ID, ES, paste0(ES_file_names[[ES]], "_", Time_step, ".tif"))
        } else if (ES_nesting[[ES]] == TRUE){
          # If ES is nesting then append time step as a dir before file path
          input_path <- file.path(input_dir, Config_ID, ES, Time_step, paste0(ES_file_names[[ES]], ".tif"))
        }
    
        # tif path structure: NCP-Time_step-Config_details$Climate_scenario.string-Config_details$Econ_scenario.string-Config_details$Pop_scenario.string-Config_details$Scenario_ID.string-full.tif’
        tif_path <- file.path(base_dir, raster_dir, paste0(ES_lower, "-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                                                Config_details$Econ_scenario.string, "-", 
                                                Config_details$Pop_scenario.string, "-", 
                                                Config_details$Scenario_ID.string, "-", names(mask), ".tif")) 
    
        png_path <- file.path(base_dir, image_dir, paste0(ES_lower, "-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                                                Config_details$Econ_scenario.string, "-", 
                                                Config_details$Pop_scenario.string, "-", 
                                                Config_details$Scenario_ID.string, "-", names(mask), ".png"))
        
        sample_path <- file.path(base_dir, sample_dir, paste0(ES_lower, "-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                                                Config_details$Econ_scenario.string, "-", 
                                                Config_details$Pop_scenario.string, "-", 
                                                Config_details$Scenario_ID.string, "-", names(mask), ".rds"))
          
        classified_area_path <- file.path(base_dir, classified_area_dir,
                            paste0(ES_lower, "-", Time_step, "-", 
                                   Config_details$Climate_scenario.string, "-", 
                                   Config_details$Econ_scenario.string, "-", 
                                   Config_details$Pop_scenario.string, "-", 
                                   Config_details$Scenario_ID.string, "-", names(mask), "-classifed_area.json"))
        

        
        # combine all paths into a named list
        time_step_paths <- list(input_path, tif_path, png_path, sample_path, classified_area_path)
        names(time_step_paths) <- c("Path", "tif_path", "png_path", "sample_path", "classified_area_path")
        return(time_step_paths)
        
        
      }))) # close loop over time_steps
    
      
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
  
  # # check if all files exist
  # if (all(ES_path_df$Exists)) {
  #   message("All ES input files exist.")
  # } else {
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
  parallel = TRUE,
  ES_layer_paths,
  ES_rescaling_dir,
  ESs_to_summarise,
  minmax_recalc = TRUE,
  report_NAs = TRUE,
  mask
  ){

  # Set parallel processing options
  if (parallel == FALSE) {
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
    future::plan(future::multisession, workers = n_workers)
  }
  
  # Create the directory to save the minmax values for each ES
  Mask_rescaling_dir <- file.path(ES_rescaling_dir, names(mask))
  if(!dir.exists(Mask_rescaling_dir)){
    dir.create(Mask_rescaling_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # loop over ESs calculating min max values for all rasters of each
  lapply(ESs_to_summarise, function(ES){

    # create an output file path for this ES
    ES_minmax_path <- file.path(Mask_rescaling_dir, paste0(tolower(ES), "_minmaxs.rds"))

    # check if file already exists & if minmax_recalc == FALSE then skip
    # otherwise if the file does not exist or minmax_recalc = TRUE then re-calculate
    if(file.exists(ES_minmax_path) & minmax_recalc == FALSE){
      cat("Global minmax file already exists for", ES, "\n")
      cat("Recalculation is set to FALSE, skipping to next ES \n")
    } else if (!file.exists(ES_minmax_path) | minmax_recalc == TRUE){

      if(!file.exists(ES_minmax_path)){
        cat("Global Minmax file does not exist for", ES, "\n")
      }
      if(minmax_recalc == TRUE){
        cat("Recalculation is set to TRUE, calculating minmax for", ES, "\n")
      }

        # Subset to ES records
        ES_layer_paths <- ES_layer_paths[ES_layer_paths$ES == ES,]

        # Inner loop over paths for ES calculating the greatest minmax values for each
        ES_minmax <- future_sapply(ES_layer_paths$Path, function(path){

          # If mask is provided, apply it to the raster
          if(names(mask) != "full_extent"){
            
            # if the mask path contains shp extension, read it as a vector
            if(grepl("\\.shp$", unlist(mask))){
              message(paste("Applying mask from shapefile:", mask))
              mask_layer <- terra::vect(x = unlist(mask))
            } else if(grepl("\\.tif$", mask)){
              message(paste("Applying mask from raster file:", mask))
              mask_layer <- terra::rast(unlist(mask))
            } else {
              stop(paste("Unsupported mask file type for:", mask))
            }
          }
          
        # Calculate min and max for current ES
        # Use try catch in case the file is missing or corrupt
        File_minmax <- tryCatch({

          # Load file
          raster <- rast(path)
          
          # If mask is provided, apply it to the raster
          if(names(mask) != "full_extent"){

            raster <- terra::mask(x = raster, 
                                  mask = mask_layer, 
                                  updatevalue = NA)
          }
          
          # Get min max
          File_minmax <- minmax(raster, compute = TRUE)

        }, error=function(e){File_minmax <- NA})

        
        cat("Finished calculating minmax for", path, "\n")
        return(File_minmax)
      }, simplify = FALSE)

      # Save minmax values for ES
      saveRDS(ES_minmax, ES_minmax_path)
    }
  })

  # if report_NAs is TRUE, check for NAs and NaNs in the minmax values
  if(report_NAs == TRUE){

    # Summarise results across ESs in a single dataframe
    All_ES_minmaxs <- lapply(ESs_to_summarise, function(ES){

      # Load the minmax values for current ES
      ES_minmax <- readRDS(file.path(Mask_rescaling_dir, paste0(tolower(ES), "_minmaxs.rds")))

      #loop over list and add details to dataframe
      ES_minmax <- lapply(ES_minmax, function(x){
        data.frame(Min = x[1], Max = x[2])
      })

      # bind list of dataframes into a single dataframe
      ES_minmax <- as.data.frame(rbindlist(ES_minmax, idcol = "Path"))
    })
    names(All_ES_minmaxs) <- ESs_to_summarise

    # bind list of dataframes into a single dataframe
    All_ES_minmaxs <- as.data.frame(rbindlist(All_ES_minmaxs, idcol = "ES"))

    # Identify NAs
    NA_indices <- which(is.na(All_ES_minmaxs$Min))

    # Get paths of the NA value
    NA_records <- ES_layer_paths[NA_indices, "Path"]

    # Check for NaN values
    NaN_indices <- which((is.nan(All_ES_minmaxs$Min)))

    # What is the path of the NaN value
    NaN_records <- ES_layer_paths[NaN_indices, "Path"]
  return(c(NA_records, NaN_records))
  }
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
    ESs_to_summarise,
    ES_rescaling_dir,
    mask
    ){
  
  # Outer loop over ESs
  ES_global_minmaxs <- lapply(ESs_to_summarise, function(ES){

    # Load the minmax values for current ES
    ES_minmax <- readRDS(file.path(ES_rescaling_dir, names(mask), paste0(tolower(ES), "_minmaxs.rds")))

    #loop over list and add details to dataframe
    ES_minmax <- lapply(ES_minmax, function(x){
      data.frame(Min = x[1], Max = x[2])
    })

    # bind list of dataframes into a single dataframe
    ES_minmax <- as.data.frame(rbindlist(ES_minmax))

    # Initialize min and max values
    ES_min <- Inf
    ES_max <- -Inf

    # Calculate global min and max for current ES
    for(i in 1:nrow(ES_minmax)){
        ES_min <- min(ES_min, ES_minmax[i, "Min"], na.rm = TRUE)
        ES_max <- max(ES_max, ES_minmax[i, "Max"], na.rm = TRUE)
    }
    return(data.frame(Min = ES_min, Max = ES_max))
  })
  names(ES_global_minmaxs) <- ESs_to_summarise

  # bind list of dataframes into a single dataframe
  ES_global_minmaxs <- as.data.frame(rbindlist(ES_global_minmaxs, idcol = "ES"))

  # Save the global minmax values
  saveRDS(ES_global_minmaxs, file = file.path(ES_rescaling_dir, names(mask), "ES_global_minmaxs.rds"))
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
save_continuous_indexed_png <- function(raster_obj, 
                                        output_path, 
                                        low_color,
                                        high_color,
                                        mid_color = NULL,
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
                                        verbose = FALSE) {
  
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
  
  # Build color palette (with optional mid_color)
  if (!is.null(mid_color)) {
    if (verbose) cat("Generating continuous palette with mid color...\n")
    color_palette <- colorRampPalette(c(low_color, mid_color, high_color))(max_colors)
  } else {
    if (verbose) cat("Generating continuous palette without mid color...\n")
    color_palette <- colorRampPalette(c(low_color, high_color))(max_colors)
  }
  
  # Create temporary PNG
  temp_png <- tempfile(fileext = ".png")
  if (verbose) cat("Creating temporary PNG...\n")
  
  png(temp_png, 
      width = width, 
      height = height, 
      res = resolution, 
      bg = background, 
      units = units)
  
  par(mar = margins)
  
  plot(raster_obj, 
       col = color_palette, 
       legend = show_legend, 
       axes = axes, 
       box = box)
  
  dev.off()
  
  # Quantize to indexed PNG
  if (verbose) cat("Converting to indexed PNG...\n")
  img <- image_read(temp_png)
  
  img_indexed <- img %>%
    image_quantize(max = max_colors, 
                   colorspace = colorspace, 
                   dither = FALSE) %>%
    image_strip()
  
  image_write(img_indexed, output_path, format = "png", depth = 8)
  
  if (verbose) cat(paste("Saved indexed PNG to:", output_path, "\n"))
  
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
normalise_layers <- function(ES_layer_paths,
                             mask,
                             ES_global_minmaxs,
                             ProjCH = ProjCH
                             ){
  
  lapply(1:nrow(ES_layer_paths), function(i){
    
    # Load the mask here because rast/vect are non-Future exportable objects
    if(names(mask) != "full_extent"){
      
      # if the mask path contains shp extension, read it as a vector
      if(grepl("\\.shp$", unlist(mask))){
        message(paste("Applying mask from shapefile:", mask))
        mask_layer <- terra::vect(x = unlist(mask))
      } else if(grepl("\\.tif$", mask)){
        message(paste("Applying mask from raster file:", mask))
        mask_layer <- terra::rast(unlist(mask))
      } else {
        stop(paste("Unsupported mask file type for:", mask))
      }
      
    }
    
    # Load the raster
    raster <- rast(ES_layer_paths[i, "Path"])
    
    # if the crs is not the same as ProjCH, reproject the raster
    if(!is.null(crs(raster)) && !isTRUE(all.equal(crs(raster), ProjCH))){
      raster <- terra::project(raster, ProjCH)
      cat("Reprojected raster to", ProjCH, "\n")
    }
    
    if(names(mask) != "full_extent"){
      # If mask is provided, apply it to the raster
      raster <- terra::mask(x = raster,
                            mask = mask_layer,
                            updatevalue = NA)
    }
    
    ES_name <- ES_layer_paths[i, "ES"]
    ES_min <- ES_global_minmaxs[1, "Min"]
    ES_max <- ES_global_minmaxs[1, "Max"]
    
    # Normalise the raster
    raster_norm <- (raster - ES_min) / (ES_max - ES_min)
    
    browser()
    
    #create dir
    dir <- dirname(ES_layer_paths[i, "tif_path"])
    if(!dir.exists(dir)){
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # save as tif using the pre-prepared path
    writeRaster(raster_norm, file = ES_layer_paths[i, "tif_path"], overwrite = TRUE)
    
    cat("rescaled layer saved to", ES_layer_paths[i, "tif_path"], "\n")
    
    # produce png img and save
    save_continuous_indexed_png(raster_obj = raster_norm,
                                output_path = ES_layer_paths[i, "png_path"],
                                low_color = "#007CDC",
                                high_color = "#FFEA2A",
                                mid_color = NULL,
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
                                verbose = FALSE)
    
    cat("rescaled layer image saved to", ES_layer_paths[i, "png_path"], "\n")
      
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
calc_summary_stats <- function(ES_layer_paths,
                             metrics = c("sum", "mean", "sd")
                             ){
  
  # Loop over sequence of ES layer paths
  ES_sum_stats <- future_lapply(1:nrow(ES_layer_paths), function(i){
    
    # use trycatch to handle missing or corrupt files
    Layer_info_stat <- tryCatch({
      
      # Load the already rescaled raster using the tif_path
      raster_norm <- rast(ES_layer_paths[i, "tif_path"])
    
      # Calculate the stats
      raster_stats <- global(raster_norm, metrics, na.rm=TRUE)
    
      # cbind raster_stats with the row of ES_layer_paths
      Layer_info_stat <- cbind(ES_layer_paths[i,], raster_stats)
      
      # return the Layer_info_stat
      return(Layer_info_stat)
    }, error=function(e){
      
      # if an error occurs, return NA for all metrics
      raster_stats <- rep(NA, length(metrics))
      names(raster_stats) <- metrics
      
      #convert to dataframe
      raster_stats <- t(data.frame(raster_stats))
      
      # cbind raster_stats with the row of ES_layer_paths
      Layer_info_stat <- cbind(ES_layer_paths[i,], raster_stats)
    })
    return(Layer_info_stat)
  })
  return(ES_sum_stats)
}



#' Calculate areas_in_classes:
#' Calculate the areas of land that fall within discrete classes of the ES provision values
#' This function calculates the areas of land that fall within discrete classes of the ES provision values
#' and saves the results to a JSON file, it is capable of calculating using different break types including
#' quantile, fisher and regular breaks that are calculated from samples taken acorss all the rescaled layers
#' called by normalise_summarize function
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param sample_size Integer: The number of samples to take from each layer for calculating global breaks
#' @param break_types Character: The types of breaks to calculate, options are "quantile", "fisher" and "regular"
#' @param save_samples Logical: Whether to save the samples to a file
Calculate_areas_in_classes <- function(ES_layer_paths,
                             sample_size = 50000,
                             save_samples = FALSE,
                             break_types = c("quantile", "fisher", "regular"),
                             save_plots = TRUE
                             ){
  
  # create an empty vector to store all samples
  all_samples <- c()
  
  # Loop over the ES_layer_paths
  for(i in 1:nrow(ES_layer_paths)){
    
    # load the rescaled raster using the tif_path
    raster_norm <- rast(ES_layer_paths[i, "tif_path"])
    
      # Take a sample for calculating the quantiles
      vals <- values(raster_norm, mat=FALSE)
      
      # remove NAs
      vals <- vals[!is.na(vals)]
      
      # vector of raster values
      samp <- sample(vals, min(sample_size, length(vals)), replace=FALSE)
    
      if(save_samples == TRUE){
    
        # save the sample
        saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
      }
      
      # add the sample to the all_samples vector
      all_samples <- c(all_samples, samp)
  }
  
  # create vector of deciles used for quantile and regular breaks
  probs <- seq(0, 1, 0.1) 
  
  # create an empty list to store breaks
  breaks <- list()
  names(breaks) <- break_types

  # loop over break_types and calculate the breaks
  for(break_type in break_types){
    
    if(break_type == "quantile"){
      # calculate quantiles
      breaks$quantile <- unique(quantile(all_samples, probs=probs, na.rm=TRUE))
    } else if(break_type == "fisher"){
      # calculate fisher natural breaks
      fisher <- classIntervals(all_samples, n = length(probs), style = "fisher")
      breaks$fisher <- fisher$brks 
    } else if(break_type == "regular"){
      # just add an empty vector for regular breaks
      breaks$regular <- rep(NA, length(probs))
    }
  }
  
  # Now loop over the layers and produce a JSON file of the % area of the map
  #that fall within each quantile
  future_lapply(1:nrow(ES_layer_paths), function(i){
    
    # Load the normalised raster using the tif_path
    raster <- rast(ES_layer_paths[i, "tif_path"])
    
    # loop over breaks and classify the raster
    for(break_type in names(breaks)){
      
      # create a raster with the breaks
      if(break_type == "quantile"){
        r_classified <- terra::classify(raster, rcl = cbind(breaks$quantile[-length(breaks$quantile)], 
                                                            breaks$quantile[-1], 1:(length(breaks$quantile)-1)))
      } else if(break_type == "fisher"){
        r_classified <- terra::classify(raster, rcl = cbind(breaks$fisher[-length(breaks$fisher)], 
                                                            breaks$fisher[-1], 1:(length(breaks$fisher)-1)))
      } else if(break_type == "regular"){
        r_classified <- terra::classify(raster, rcl = cbind(probs[-length(probs)], probs[-1], 1:(length(probs)-1)))
      }
      
      # get frequency of each bin (counts of cells)
      freq_tbl <- freq(r_classified)  # terra version of raster::freq
      freq_tbl <- as.data.frame(freq_tbl)
      
      # total non-NA cells
      total_cells <- sum(freq_tbl$count)
      
      # percentage per bin
      freq_tbl$percent <- 100 * freq_tbl$count / total_cells
      
      # create labels for the bins
      make_labels <- function(breaks, n_bins, min_value = 0) {
        lower <- c(min_value, breaks[2:n_bins])
        upper <- breaks[2:(n_bins+1)]
        paste0(round(lower, 2), "–", round(upper, 2))
      }
      
      # Assign labels to the frequency table
      freq_tbl$bin_label <- make_labels(breaks[[break_type]], n_bins = length(breaks[[break_type]]) - 1)
      
      # add method label
      freq_tbl$method <- break_type
      
      # save the frequency table as a JSON file
      json_data <- toJSON(setNames(as.list(freq_tbl$bin_label), freq_tbl$percent), pretty = TRUE)
      
      # modify the classified area path to include the break type
      ES_layer_paths[i, "classified_area_path"] <- gsub(".json", paste0("_", break_type, ".json"), 
                                                        ES_layer_paths[i, "classified_area_path"])
      
      # save the json data to the path_classified_area
      write(json_data, file = ES_layer_paths[i, "classified_area_path"])
      
      # if save_plots is TRUE, save a bar chart of the frequency table
      if(save_plots == TRUE){
        
        # create a bar chart of the frequency table
        bar_plot <- ggplot(freq_tbl, aes(x = bin_label, y = percent)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          labs(title = paste("Area in Classes for", ES_layer_paths[i, "ES"], "using", break_type, "breaks"),
               x = "Class",
               y = "Percentage of Area") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        # save the plot as a PNG file
        ggsave(filename = gsub(".json", paste0("_", break_type, ".png"), 
                               ES_layer_paths[i, "classified_area_path"]),
               plot = bar_plot,
               width = 10, height = 6)
      }
      
    }
  })
  cat("Finished calculating breaks for all layers of current ES and saving classified area JSON files.\n")
}


#' normalise_and_summarize:
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
#' @param Recalc_rescaled_layers Logical: If TRUE, the function will recalculate the rescaled layers
#' If FALSE, the function will skip the recalculation if the rescaled layers already exist
#' @param Recalc_summary Logical: If TRUE, the function will recalculate the summary statistics
#' If FALSE, the function will skip the recalculation if the summary statistics already exist
#' @param Recalc_classified_area Logical: If TRUE, the function will recalculate the classified area
#' If FALSE, the function will skip the recalculation if the classified area already exists
#' @param Sim_ctrl_tbl_path Character: The path to the simulation control table
#' @param Summarize_across_ES Logical: If TRUE, the function will summarize the results across all ESs
#' @param Rescale_cross_ES_results Logical: If TRUE, the function will rescale the summary results across all ESs
#' @param Parallel Logical: If TRUE, the function will run in parallel using the number of workers specified in the config file
#' @param mask Character: The path to the mask to apply to the layers
#' @param sample_size Integer: The number of samples to take from each layer for calculating global breaks
#' @param classified_area_break_types Character: The types of breaks to calculate for classified area, options are "quantile", "fisher" and "regular"
#' @param save_classified_area_plots Logical: Whether to save the classified area plots
#' @param save_break_samples Logical: Whether to save the samples used for calculating breaks
#' @return: NULL
#' 
#' @export
#
normalise_and_summarize <- function(
    ProjCH = "+proj=somerc +init=epsg:2056",
    metrics = c("sum", "mean", "sd"),
    ES_layer_paths,
    ESs_to_summarise,
    ES_rescaling_dir,
    ES_summary_stats_dir,
    ES_summarisation_dir,
    Recalc_rescaled_layers = FALSE,
    Recalc_summary = TRUE,
    Recalc_classified_area = FALSE,
    Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
    Summarize_across_ES = TRUE,
    Rescale_cross_ES_results = TRUE,
    Parallel = TRUE,
    mask,
    sample_size = 50000,
    classified_area_break_types = c("quantile", "fisher", "regular"),
    save_classified_area_plots = TRUE,
    save_break_samples = FALSE
  ){

  # Set parallel processing options
  if (Parallel == FALSE) {
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
    future::plan(future::multisession, workers = n_workers)
  }
  
  # Load back in the global minmax values
  ES_global_minmaxs <- readRDS(file.path(ES_rescaling_dir, names(mask), "ES_global_minmaxs.rds"))
  
  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)
  
  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI_")
  
  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "GREX", "GR_EX")

  Mask_summary_stats_dir <- file.path(ES_summary_stats_dir, names(mask))
  
  # if the directory does not exist, create it
  if(!dir.exists(Mask_summary_stats_dir)){
    dir.create(Mask_summary_stats_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Loop over the ESs specified and calculate the summary metrics specified for the
  # rescaled values for each time step under all configurations
  # for each configuration and time-step
  lapply(ESs_to_summarise, function(ES){
    
    # convert ES to lowercase
    ES <- tolower(ES)
    
    # subset to ES records
    ES_layer_paths <- ES_layer_paths[ES_layer_paths$ES == ES,]
    
    # subset the global minmaxs to the current ES
    ES_global_minmaxs <- ES_global_minmaxs[ES_global_minmaxs$ES == ES,]
    
    # 1. Rescale the ES layers
    if(Recalc_rescaled_layers == TRUE){
      
      # Apply the normalisation function to the ES layers
      normalise_layers(ES_layer_paths = ES_layer_paths,
                       mask = mask,
                      ES_global_minmaxs = ES_global_minmaxs)
      
    } else if(Recalc_rescaled_layers == FALSE){
      cat("Recalculation of rescaled layers is set to FALSE, identifying which layers already exist for", ES, "\n")
      
      # Count how many of the rescaled layers for this ES exist and how many
      # are missing
      ES_layer_paths$Exists <- sapply(ES_layer_paths$tif_path, file.exists)
      
      cat("Number of rescaled raster layers for", ES, "that exist:", sum(ES_layer_paths$Exists), "\n")
      cat("Number of rescaled raster layers for", ES, "that are missing:", sum(!ES_layer_paths$Exists), "\n")
      
      # Remove the existing layers from the dataframe
      ES_layer_paths_remaining <- ES_layer_paths[!ES_layer_paths$Exists,]
      
      cat("Existing layers have been removed from the list to be rescaled, performing normalization of remaining layers \n")
      
      # If there are no layers to rescale, skip to the next ES
      if(nrow(ES_layer_paths_remaining) == 0){
        cat("Recalculation of rescaled layers is set to FALSE and all layers already exist for", ES, "skipping to next step \n")
      } else if (nrow(ES_layer_paths_remaining) > 0){
        
        cat("Recalculation of rescaled layers is set to FALSE, normalising remaining layers for", ES, "\n")
        
        # If there are layers to rescale, apply the normalisation function
        normalise_layers(ES_layer_paths = ES_layer_paths_remaining,
                         mask = mask,
                         ES_global_minmaxs = ES_global_minmaxs)
      }
    
    }
    
    # 2. Calculate summary statistics for the rescaled layers
    # create path to save the combined results for this ES
    ES_sum_stats_path <- file.path(Mask_summary_stats_dir, paste0(tolower(ES), "_summary_stats.rds"))
    
    # check if the summary file already exists & if Recalc_summary == FALSE then skip
    # otherwise if the file does not exist or Recalc_summary = TRUE then re-calculate
    if(file.exists(ES_sum_stats_path) & Recalc_summary == FALSE){
      cat("Summary stats file already exists for", ES, "\n")
      cat("Recalculation of summary is set to FALSE, skipping to next step \n")
    } else if (!file.exists(ES_sum_stats_path) | Recalc_summary == TRUE){
      
      if(!file.exists(ES_sum_stats_path)){
        cat("Summary stats file does not exist for", ES, "\n")
      }
      if(Recalc_summary == TRUE){
        cat("Recalculation of summary is set to TRUE, calculating Summary stats for", ES, "\n")
      }
      
      # Apply function to calculate summary statistics for the rescaled layers
      ES_sum_stats <- calc_summary_stats(ES_layer_paths,
                             metrics = metrics
                             )
      # Save the results for the current ES
      saveRDS(ES_sum_stats, file = ES_sum_stats_path)
      cat("Summary stats for", ES, "have been saved to", ES_sum_stats_path, "\n")
    }
      
    browser()
    # 3. Classify the rescaled layers and calculate the area of each class
    if(Recalc_classified_area == TRUE){
      
      # Calculate the areas in classes for the rescaled layers
      Calculate_areas_in_classes(ES_layer_paths = ES_layer_paths,
                                 sample_size = sample_size,
                                 save_samples = save_break_samples,
                                 break_types = classified_area_break_types,
                                 save_plots = save_classified_area_plots)
      cat("Areas in classes for", ES, "have been calculated and saved to JSON files \n")
    }
  }) # close loop over ESs_to_summarise

  # if parallel == TRUE set back to sequential processing
  if(Parallel == TRUE){
    future::plan(future::sequential)
  }
  
  # 4. Summarise the results across all ESs
  if(Summarize_across_ES == TRUE){
    
    # Produce a single df summarising the results across all ESs
    ES_sum_stats <- lapply(ESs_to_summarise, function(ES){
      
      # Load the summary stats for the current ES
      ES_sum_stats <- readRDS(file.path(Mask_summary_stats_dir, paste0(tolower(ES), "_summary_stats.rds")))
      
      # Bind the list of dataframes for each time point together
      ES_sum_stats <- do.call(rbind, ES_sum_stats)
      
      # remove the exists column
      ES_sum_stats$Exists <- NULL
      
      return(ES_sum_stats)
    })
    
    #rbind list of dataframes into a single dataframe
    ES_sum_stats <- do.call(rbind, ES_sum_stats)
    
    # Use the simulation control table to add the scenario_ID to the ES_sum_stats
    ES_sum_stats$Scenario <- sapply(ES_sum_stats$Config_ID, function(x) Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == x, "Scenario_ID.string"])
    
    # Convert metric columns to numeric
    ES_sum_stats[,metrics] <- lapply(ES_sum_stats[,metrics], as.numeric)
    
    # create dir for saving summarisation results
    Mask_summarisation_dir <- file.path(ES_summarisation_dir, names(mask))
    if(!dir.exists(Mask_summarisation_dir)){
      dir.create(Mask_summarisation_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save the ES_sum_stats
    saveRDS(ES_sum_stats, file = file.path(Mask_summarisation_dir, "ES_summary_stats.rds"))
    
    cat("Summary stats across all ES have been calculated and saved \n")
    
    # if Rescale_results is TRUE, then normalise the results and save a copy
    if(Rescale_cross_ES_results == TRUE){
      cat("rescaled_summary_stats is set to TRUE, normalising the summary stats \n")
      
      # Normalise the summary stats
      ES_sum_stats_norm <- lapply(ESs_to_summarise, function(ES){
        
        # Get summary values for current ES
        ES_sum_stats_ES <- ES_sum_stats[ES_sum_stats$ES == ES,]
        
        # Loop over the metrics
        for(metric in metrics){
          
          # Initialize min and max values
          ES_min <- min(ES_sum_stats_ES[[metric]], na.rm = TRUE)
          ES_max <- max(ES_sum_stats_ES[[metric]], na.rm = TRUE)
          
          # Rescale the summary values
          ES_sum_stats_ES[[paste0(metric, "_rescale")]] <- (ES_sum_stats_ES[[metric]] - ES_min) / (ES_max - ES_min)
          
          # Remove the original summary values
          ES_sum_stats_ES[[metric]] <- NULL
          
          # Rename the rescaled summary values
          names(ES_sum_stats_ES)[names(ES_sum_stats_ES) == paste0(metric, "_rescale")] <- metric
        }
        
        return(ES_sum_stats_ES)
      })
      names(ES_sum_stats_norm) <- ESs_to_summarise
      
      # bind the list of dataframes into a single dataframe
      ES_sum_stats_norm <- as.data.frame(rbindlist(ES_sum_stats_norm))
      
      # Save the rescaled summary stats
      saveRDS(ES_sum_stats_norm, file = file.path(Mask_summarisation_dir, "ES_summary_stats_rescaled.rds"))
    }
  }
  
  cat("Finished normalising and summarising the ES layers \n")
}


#' Summarize_for_masks
#' Upper level function to summarise the ESs for a list of supplied masks
#' #' @param config List: Configuration list containing the ESs to summarise,
#' #' ES nesting, ES file names, and other parameters
#' #' @param web_platform_dir Character: The directory where the web platform files are stored
#' #' @param ProjCH string of CRS
#' #' @param ES_rescaling_dir Character: The directory where the rescaled ES layers are stored
#' #' @param masks List: A named list of masks to apply to the ES layers if the
#' #  name is full_extent then the list item can be empty and the summary will be performed on the entire map


masks <- list("canton" = file.path(Mask_dir, "Canton_mask.shp"))
i=1


Summarise_for_masks <- function(
    config = config,
    web_platform_dir,
    ProjCH,
    ES_rescaling_dir = ES_rescaling_dir,
    ES_summary_stats_dir = ES_summary_stats_dir,
    ES_summarisation_dir = ES_summarisation_dir,
    masks = list("canton" = file.path(Mask_dir, "Canton_mask.shp"))
){
  
  
  # loop over the masks
  for(i in seq_along(masks)){
    
    cat(paste0("Processing ES outputs for mask: ", names(masks)[i], "\n"))
    
    # apply function to prepare DF
    ES_layer_paths <- prepare_analysis_df(
      input_dir = "F:/KB-outputs/ncp_output",
      image_dir = "map_images",
      raster_dir = "raster_data",
      classified_area_dir = "chart_data/classified_area",
      base_dir = web_platform_dir,
      Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
      ProjCH = ProjCH,
      ESs_to_summarise = config$ESs_to_summarise,
      ES_nesting = config$ES_nesting,
      ES_file_names = config$ES_file_names,
      mask = masks[i]
    )
    
    # Calculate minmaxs for this mask
    # calc_minmaxs(
    #   parallel = config$Parallel,
    #   ES_layer_paths = ES_layer_paths,
    #   ES_rescaling_dir = ES_rescaling_dir,
    #   ESs_to_summarise = config$ESs_to_summarise,
    #   minmax_recalc = TRUE,
    #   report_NAs = TRUE,
    #   mask = masks[i]
    # )
    # 
    # # Calculate global minmaxs for this mask
    # calc_global_minmaxs(
    #   ESs_to_summarise = config$ESs_to_summarise,
    #   ES_rescaling_dir = ES_rescaling_dir,
    #   mask = masks[i]
    # )
    
    # Normalise and summarise the ESs for this mask
    normalise_and_summarize(
      ProjCH = ProjCH,
      metrics = c("sum", "mean", "sd"),
      ES_layer_paths = ES_layer_paths,
      ESs_to_summarise = config$ESs_to_summarise,
      ES_rescaling_dir = ES_rescaling_dir,
      ES_summary_stats_dir = ES_summary_stats_dir,
      ES_summarisation_dir = ES_summarisation_dir,
      Recalc_rescaled_layers = TRUE,
      Recalc_summary = TRUE,
      Recalc_classified_area = FALSE,
      Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
      Summarize_across_ES = TRUE,
      Rescale_cross_ES_results = TRUE,
      Parallel = TRUE,
      mask = masks[i],
      sample_size = 50000,
      classified_area_break_types = c("quantile", "fisher", "regular"),
      save_classified_area_plots = TRUE,
      save_break_samples = FALSE
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
  config$InputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR,
                               bash_vars$ES_OUTPUT_BASE_DIR)
}
# if OutputDir is not set, use bash variable SUMMARISATION_OUTPUT_DIR
if (is.null(config$OutputDir) || config$OutputDir == "") {
  config$OutputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR,
                                bash_vars$SUMMARISATION_OUTPUT_DIR)
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
if (!dir.exists(config$OutputDir)) {
  dir.create(config$OutputDir, recursive = TRUE)
}

cat("Performing summarisation of ES outputs \n")

cat("Working directory is:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")

### Start by creating sub directories

# dir for saving results
web_platform_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform"
# if dir doesn't exist, create it
if(!dir.exists(web_platform_dir)){
  dir.create(web_platform_dir, recursive = TRUE, showWarnings = FALSE)
}

# Mask directory
Mask_dir <- file.path(config$OutputDir, "Masks")
# if dir doesn't exist, create it
if(!dir.exists(Mask_dir)){
  dir.create(Mask_dir, recursive = TRUE, showWarnings = FALSE)
}

# Dir for the processed ES layers 
# Sub-dir for the summarisation results
ES_summarisation_dir <- file.path(config$OutputDir, "ES_summarisation")
if (!dir.exists(ES_summarisation_dir)) {
  dir.create(ES_summarisation_dir, recursive = TRUE)
}

# Sub-dir for the normalisation results
ES_rescaling_dir <- file.path(ES_summarisation_dir, "ES_rescaling")
if (!dir.exists(ES_rescaling_dir)) {
  dir.create(ES_rescaling_dir, recursive = TRUE)
}

# Sub-dir for the summary_stats results
ES_summary_stats_dir <- file.path(ES_summarisation_dir, "ES_summary_stats")
if (!dir.exists(ES_summary_stats_dir)) {
  dir.create(ES_summary_stats_dir, recursive = TRUE)
}

# Dir to store results of change analysis
ES_change_dir <- file.path(ES_summarisation_dir, "ES_change")
if (!dir.exists(ES_change_dir)) {
  dir.create(ES_change_dir, recursive = TRUE)
}

# Sub-dir for Change_summary_stats and 
ES_change_summary_stats_dir <- file.path(ES_change_dir, "summary_stats")
if (!dir.exists(ES_change_summary_stats_dir)) {
  dir.create(ES_change_summary_stats_dir, recursive = TRUE)
}

# Sub-dir for rasters summing change over time
ES_change_sum_rasters_dir <- file.path(ES_change_dir, "sum_change_rasters")
if (!dir.exists(ES_change_sum_rasters_dir)) {
  dir.create(ES_change_sum_rasters_dir, recursive = TRUE)
}

# Sub-dir to store results of SSIM pattern analysis
ES_SSIM_dir <- file.path(ES_summarisation_dir, "ES_SSIM")
if (!dir.exists(ES_SSIM_dir)) {
  dir.create(ES_SSIM_dir, recursive = TRUE)
}


