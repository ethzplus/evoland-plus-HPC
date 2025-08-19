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
           "tidyterra", "jsonlite", "magick", "grDevices")
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

  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)
  
  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI-")
  
  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "GREX", "GR-EX")
  
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
        tif_path <- file.path(base_dir, raster_dir, paste0(ES, "-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                                                Config_details$Econ_scenario.string, "-", 
                                                Config_details$Pop_scenario.string, "-", 
                                                Config_details$Scenario_ID.string, "-", names(mask), ".tif")) 
    
        png_path <- file.path(base_dir, image_dir, paste0(ES, "-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                                                Config_details$Econ_scenario.string, "-", 
                                                Config_details$Pop_scenario.string, "-", 
                                                Config_details$Scenario_ID.string, "-", names(mask), ".png"))
          
        classified_area_path <- file.path(base_dir, classified_area_dir,
                            paste0(ES, "-", Time_step, "-", 
                                   Config_details$Climate_scenario.string, "-", 
                                   Config_details$Econ_scenario.string, "-", 
                                   Config_details$Pop_scenario.string, "-", 
                                   Config_details$Scenario_ID.string, "-", names(mask), "-classifed_area.json"))
        

        
        # combine all paths into a named list
        time_step_paths <- list(input_path, tif_path, png_path, classified_area_path)
        names(time_step_paths) <- c("Path", "es_path_tif", "es_path_png", "es_path_classifed_area")
        return(time_step_paths)
        
        
      }))) # close loop over time_steps
    
      
      #Add column for ES
      Config_time_paths$ES <- ES
      
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
    ES_minmax_path <- file.path(Mask_rescaling_dir, paste0(ES, "_minmaxs.rds"))

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
          if(names(mask) != "Full_extent"){
            
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
          if(names(mask) != "Full_extent"){

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
      ES_minmax <- readRDS(file.path(Mask_rescaling_dir, paste0(ES, "_minmaxs.rds")))

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
    ES_minmax <- readRDS(file.path(ES_rescaling_dir, names(mask), paste0(ES, "_minmaxs.rds")))

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


#' normalise_and_summarize:
#' 
#' Calculate the sum, mean and standard deviation of the normalised values
#'  for each time step of each configuration for the ESs specificed 
#' 
#' This function calculates the sum, mean and standard deviation of the normalised values
#' for each time step of each configuration for the ESs specificed and saves the results
#' to an RDS file
#' 
#' @param metrics Vector of metrics to calculate must use built-in function arguments from terra::global
#' i.e. "max", "min", "mean", "sum", "range", "rms" (root mean square), "sd",
#'  "std" (population sd, using n rather than n-1),
#' @param ES_layer_paths Dataframe: containing info for each ES layer
#' i.e. config_ID, time step, raw layer path and path to save rescaled layer
#' @param ESs_to_summarise Vector of ESs to summarise
#' @param ES_rescaling_dir Directory to load the global minmaxs values used for normalization 
#' @param ES_summary_stats_dir Directory to save the summary statistics for individual ESs
#' @param ES_summarisation_dir Directory to save the summary statistics for all ESs
#' @param Recalc_summary Logical indicating whether to recalculate the summary statistics for each ES
#' @param Recalc_rescaled_layers Logical indicating whether to recalculate the normalised layers
#' @param Save_rescaled_layers Logical indicating whether to save the normalised layers
#' @param Sim_ctrl_tbl Dataframe: containing the simulation information
#' @param Rescale_results Logical indicating whether to normalise the results
#' @param Parallel Logical indicating whether to use parallel processing
#' 
#' @export
#
normalise_and_summarize <- function(
    metrics = c("sum", "mean", "sd"),
    ES_layer_paths,
    ESs_to_summarise,
    ES_rescaling_dir,
    ES_summary_stats_dir,
    ES_summarisation_dir,
    Recalc_summary = TRUE,
    Recalc_rescaled_layers = FALSE,
    Save_rescaled_layers = FALSE,
    Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
    Rescale_results = TRUE,
    Parallel = TRUE,
    mask,
    sample_size = 50000
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
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI-")
  
  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "GREX", "GR-EX")

  # Loop over the ESs specified and calculate the summary metrics specified for the
  # rescaled values for each time step under all configurations
  # for each configuration and time-step
  lapply(ESs_to_summarise, function(ES){
    
    # subset to ES records
    ES_layer_paths <- ES_layer_paths[ES_layer_paths$ES == ES,]
    
    Mask_summary_stats_dir <- file.path(ES_summary_stats_dir, names(mask))
    # if the directory does not exist, create it
    if(!dir.exists(Mask_summary_stats_dir)){
      dir.create(Mask_summary_stats_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # create path to save the combined results for this ES
    ES_sum_stats_path <- file.path(Mask_summary_stats_dir, paste0(ES, "_summary_stats.rds"))
    
    # check if the summary file already exists & if Recalc_summary == FALSE then skip
    # otherwise if the file does not exist or Recalc_summary = TRUE then re-calculate
    if(file.exists(ES_sum_stats_path) & Recalc_summary == FALSE){
      cat("Summary stats file already exists for", ES, "\n")
      cat("Recalculation of summary is set to FALSE, skipping to next ES \n")
    } else if (!file.exists(ES_sum_stats_path) | Recalc_summary == TRUE){
      
      if(!file.exists(ES_sum_stats_path)){
        cat("Summary stats file does not exist for", ES, "\n")
      }
      if(Recalc_summary == TRUE){
        cat("Recalculation of summary is set to TRUE, calculating Summary stats for", ES, "\n")
      }
      
      #if Recalc_rescaled_layers is FALSE, then identify which layers exist
      #already and remove them from the df
      if(Recalc_rescaled_layers == FALSE){
        
        cat("Recalculation of rescaled layers is set to FALSE, identifying which layers already exist for", ES, "\n")
        
        # Count how many of the rescaled layers for this ES exist and how many
        # are missing
        ES_layer_paths$Exists <- sapply(ES_layer_paths$es_path_tif, file.exists)
        
        cat("Number of rescaled raster layers for", ES, "that exist:", sum(ES_layer_paths$Exists), "\n")
        cat("Number of rescaled raster layers for", ES, "that are missing:", sum(!ES_layer_paths$Exists), "\n")
        
        # Remove the existing layers from the dataframe
        ES_layer_paths <- ES_layer_paths[!ES_layer_paths$Exists,]
        
        cat("Existing layers have been removed from the list to be rescaled, performing normalization of remaining layers \n")
      } else if(Recalc_rescaled_layers == TRUE){
        cat("Recalculation of rescaled layers is set to TRUE, normalising all layers in ES_layer_paths \n")
      }
      
      if(nrow(ES_layer_paths) == 0){
        cat("Recalculation of rescaled layers is set to FALSE and all layers already exist for", ES, "skipping to next ES \n")
      }else{
        
        # Loop over sequence of ES layer paths, normalise the layer values 
        # with the global min/max values and then calculate sum, mean, median, mode
        # std of the rescaled values. 
        ES_sum_stats <- future_lapply(1:nrow(ES_layer_paths), function(i){
          
          # Load the mask here because rast/vect are non-Future exportable objects
          if(names(mask) != "Full_extent"){
            
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
          
          # use trycatch to handle missing or corrupt files
          Layer_info_stat <- tryCatch({
            
            # Load the raster
            raster <- rast(ES_layer_paths[i, "Path"])
            
            if(names(mask) != "Full_extent"){
              # If mask is provided, apply it to the raster
              raster <- terra::mask(x = raster, 
                                    mask = mask_layer, 
                                    updatevalue = NA)
            }
            
            ES_name <- ES_layer_paths[i, "ES"]
            ES_min <- ES_global_minmaxs[ES_global_minmaxs$ES == ES, "Min"]
            ES_max <- ES_global_minmaxs[ES_global_minmaxs$ES == ES, "Max"]
            
            # Normalise the raster
            raster_norm <- (raster - ES_min) / (ES_max - ES_min)
            
            # if Save_rescaled_layers is TRUE, save the rescaled raster
            if(Save_rescaled_layers == TRUE){
              
              cat("Save_rescaled_layers is set to TRUE \n")
              
              #create dir
              dir <- dirname(ES_layer_paths[i, "es_path_tif"])
              if(!dir.exists(dir)){
                dir.create(dir, recursive = TRUE, showWarnings = FALSE)
              }
              
              # save as tif using the pre-prepared path
              writeRaster(raster_norm, file = ES_layer_paths[i, "es_path_tif"], overwrite = TRUE)
              
              cat("rescaled layer saved to", ES_layer_paths[i, "es_path_tif"], "\n")
              
              # produce png img and save
              save_continuous_indexed_png(raster_obj = raster_norm, 
                                          output_path = ES_layer_paths[i, "es_path_png"],
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
              
              cat("rescaled layer image saved to", ES_layer_paths[i, "es_path_png"], "\n")
            }
            
            # Calculate the stats
            raster_stats <- global(raster_norm, metrics, na.rm=TRUE)
            
            # Take a sample for calculating the quantiles
            vals <- values(raster_norm, mat=FALSE)
            
            # remove NAs
            vals <- vals[!is.na(vals)]
            
            # vector of raster values
            samp <- sample(vals, min(sample_size, length(vals)), replace=FALSE) 
            
            # save the sample as an RDS using the mask name, ES and Config_ID
            sample_path <- file.path(ES_rescaling_dir, names(mask), "quantile_samples",
                                     paste0(ES, "_sample_", ES_layer_paths[i, "Config_ID"], ".rds"))
            
            # create the directory if it does not exist
            sample_dir <- dirname(sample_path)
            if(!dir.exists(sample_dir)){
              dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
            }
            
            # save the sample
            saveRDS(samp, file = sample_path)
            
            # cbind raster_stats with the row of ES_layer_paths
            Layer_info_stat <- cbind(ES_layer_paths[i,], raster_stats)
            
            # add the sample path to the Layer_info_stat
            Layer_info_stat$Sample_path <- sample_path
            
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

            # add the sample path as NA
            Layer_info_stat$Sample_path <- NA
          })
          return(Layer_info_stat)
        })
        
        # Save the results for the current ES
        saveRDS(ES_sum_stats, file = ES_sum_stats_path)
        cat("Summary stats for", ES, "have been saved to", ES_sum_stats_path, "\n")
      
        # loop over the non-NA Sample_paths and calculate the quantiles
        ES_sample_paths <- ES_sum_stats[!is.na(ES_sum_stats)]
        
        all_samples <- c()
        for(sample_path in ES_sample_paths){
          
          # Load the sample
          sample <- readRDS(sample_path)
          
          # add the sample to the all_samples vector
          all_samples <- c(all_samples, sample)
        }
        
        # Calculate the quantiles for the combined sample
        probs <- seq(0, 1, 0.1)  # deciles, change as needed
        global_qtiles <- quantile(all_samples, probs=probs, na.rm=TRUE)
        
        # Save the global quantiles to a file
        global_qtiles_path <- file.path(Mask_summary_stats_dir, paste0(ES, "_global_quantiles.rds"))
        saveRDS(global_qtiles, global_qtiles_path)
        
        cat("Saving global quantiles for", ES, "to", global_qtiles_path, "\n")
      
        browser()
        # Now loop over the layers and produce a JSON file of the % area of the map
        #that fall within each quantile
        lapply(1:nrow(ES_layer_paths), function(i){

          # Load the raster
          raster <- rast(ES_layer_paths[i, "Path"])

          # If mask is provided, apply it to the raster
          if(names(mask) != "Full_extent"){
            raster <- terra::mask(x = raster, 
                                  mask = mask_layer, 
                                  updatevalue = NA)
          }

          # Calculate the area of each quantile as a % of the area of NON-NA cells in the raster
          Non_na_area <- sum(!is.na(values(raster)), na.rm = TRUE)
          quantile_areas <- sapply(global_qtiles, function(qtile){
            sum(values(raster) <= qtile, na.rm = TRUE) / Non_na_area * 100
          })

          # Create a data frame with the quantile areas
          quantile_df <- data.frame(Quantile = names(global_qtiles), Area = quantile_areas)
          
          # convert the df to json
          json_data <- toJSON(setNames(as.list(quantile_df$Quantile), quantile_df$Area), pretty = TRUE)
          
          # save the json data to the mask_path_data
          write(json_data, file = ES_layer_paths[i, "es_path_classified_area"])
          cat("Quantile areas for", ES_layer_paths[i, "Path"], "have been saved to", ES_layer_paths[i, "es_path_classified_area"], "\n")
          browser()
        })
        
        
        } # close if statement for nrow(ES_layer_paths) == 0
      }# close if statement for Recalc_rescaled_layers
  }) # close loop over ESs_to_summarise

  # if parallel == TRUE set back to sequential processing
  if(Parallel == TRUE){
    future::plan(future::sequential)
  }
  

  
  
  
  

  # Produce a single df summarising the results across all ESs
  ES_sum_stats <- lapply(ESs_to_summarise, function(ES){

    # Load the summary stats for the current ES
    ES_sum_stats <- readRDS(file.path(Mask_summary_stats_dir, paste0(ES, "_summary_stats.rds")))

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

  cat("Summary stats have been calculated and saved \n")

  # if Rescale_results is TRUE, then normalise the results and save a copy
  if(Rescale_results == TRUE){
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

#' save_continuous_indexed_png
#' # Save a continuous indexed PNG from a raster object with specified colors
#' #' @param raster_obj A RasterLayer or RasterStack object to be saved as PNG
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



#' Summarize_for_masks
#' Upper level function to summarise the ESs for a list of supplied masks
#' #' @param config List: Configuration list containing the ESs to summarise,
#' #' ES nesting, ES file names, and other parameters
#' #' @param web_platform_dir Character: The directory where the web platform files are stored
#' #' @param ProjCH string of CRS
#' #' @param ES_rescaling_dir Character: The directory where the rescaled ES layers are stored
#' #' @param masks List: A named list of masks to apply to the ES layers if the
#' #  name is Full_extent then the list item can be empty and the summary will be performed on the entire map


masks <- list("full" = file.path(Mask_dir, "Canton_mask.shp"))
i=1


Summarise_for_masks <- function(
    config = config,
    web_platform_dir,
    ProjCH,
    ES_rescaling_dir = ES_rescaling_dir,
    ES_summary_stats_dir = ES_summary_stats_dir,
    ES_summarisation_dir = ES_summarisation_dir,
    masks = list("full" = file.path(Mask_dir, "Canton_mask.shp"))
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
    calc_minmaxs(
      parallel = config$Parallel,
      ES_layer_paths = ES_layer_paths,
      ES_rescaling_dir = ES_rescaling_dir,
      ESs_to_summarise = config$ESs_to_summarise,
      minmax_recalc = TRUE,
      report_NAs = TRUE,
      mask = masks[i]
    )
    
    # Calculate global minmaxs for this mask
    calc_global_minmaxs(
      ESs_to_summarise = config$ESs_to_summarise,
      ES_rescaling_dir = ES_rescaling_dir,
      mask = masks[i]
    )
    
    # Normalise and summarise the ESs for this mask
    normalise_and_summarize(
      metrics = c("sum", "mean", "sd"),
      ES_layer_paths = ES_layer_paths,
      ESs_to_summarise = config$ESs_to_summarise,
      ES_rescaling_dir = ES_rescaling_dir,
      ES_summary_stats_dir = ES_summary_stats_dir,
      ES_summarisation_dir = ES_summarisation_dir,
      Recalc_summary = TRUE,
      Recalc_rescaled_layers = TRUE,
      Save_rescaled_layers = TRUE,
      Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
      Rescale_results = TRUE,
      Parallel = TRUE,
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


