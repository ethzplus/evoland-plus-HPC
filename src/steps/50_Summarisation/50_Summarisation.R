#' Calculate various summary metrics across the configurations
#'
#' Input: Folder of outputs from NCP simulations:
#'    Normalise the NCP outputs and calculate a series of different summary
#'    metrics from them, save results as rds files to be used in visualisations
#'    and subsequent analysis
#' Output: Summary measures of NCP provision across configurations. 
#'
#' @environment Summarisation
#' @config $FUTURE_EI_CONFIG_FILE (yaml file)
#' @date 2024-10-15
#' @author Benjamin Black
#'
#'
#' @docType script
#'

# # FOR TESTING ONLY
# # Create an environment variable for the config file
# Sys.setenv(FUTURE_EI_CONFIG_FILE = "config.yml")

# Load libraries
packs <- c("stringr", "terra", "future", "future.apply",
           "data.table", "tidyr", "yaml", "dplyr")
invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

## Functions

#' calc_minmaxs:
#'
#' Calculate the minmax values for all configurations for each NCP
#'
#' This function calculates the minimum and maximum values in the rasters 
#' for all time points for all configurations for each NCP  and saves 
#' the results for each NCP to an RDS file
#'
#'
#' @param parallel: logical: If TRUE, the function will run in parallel
#' using the number of workers specified in the config file. If FALSE, the
#' function will run sequentially
#' @param NCP_layer_paths: data.frame: A df containing info for each NCP layer
#' i.e. config_ID, time step, raw layer path and path to save normalized layer
#' @param NCP_normalisation_dir: character: The directory to save the minmax
#' values for each NCP
#' @param NCPs_to_summarise: character: A vector of NCPs to calculate minmax
#' values for
#' @param minmax_recalc: logical: If TRUE, the function will recalculate the
#' minmax values for each NCP even if the file already exists. If FALSE, the
#' function will skip the calculation if the file already exists
#' @param report_NAs: logical: If TRUE, the function will check for NAs and NaNs
#' in the minmax values. If any are found, the function will return a vector of
#' the NCP paths that have NAs or NaNs as the min or max values indicating
#' there may be a problem with the layers
#'
#'@return: vector: If report_NAs is TRUE, the function will return a vector of
#' the NCP paths that have NAs or NaNs as the min or max values indicating 
#' there may be a problem with the layers
#' @export

calc_minmaxs <- function(
  parallel = TRUE,
  NCP_layer_paths,
  NCP_normalisation_dir,
  NCPs_to_summarise,
  minmax_recalc = TRUE,
  report_NAs = TRUE  
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

  # loop over NCPs calculating min max values for all rasters of each
  lapply(NCPs_to_summarise, function(NCP){
  
    # create an output file path for this NCP
    NCP_minmax_path <- file.path(NCP_normalisation_dir, paste0(NCP, "_minmaxs.rds"))
  
    # check if file already exists & if minmax_recalc == FALSE then skip
    # otherwise if the file does not exist or minmax_recalc = TRUE then re-calculate
    if(file.exists(NCP_minmax_path) & minmax_recalc == FALSE){
      cat("Global minmax file already exists for", NCP, "\n")
      cat("Recalculation is set to FALSE, skipping to next NCP \n")
    } else if (!file.exists(NCP_minmax_path) | minmax_recalc == TRUE){
    
      if(!file.exists(NCP_minmax_path)){
        cat("Global Minmax file does not exist for", NCP, "\n")
      }
      if(minmax_recalc == TRUE){
        cat("Recalculation is set to TRUE, calculating minmax for", NCP, "\n")
      }
    
        # Subset to NCP records
        NCP_layer_paths <- NCP_layer_paths[NCP_layer_paths$NCP == NCP,]
    
        # Inner loop over paths for NCP calculating the greatest minmax values for each
        NCP_minmax <- future_sapply(NCP_layer_paths$Path, function(path){
      
        # Calculate min and max for current NCP
        # Use try catch in case the file is missing or corrupt
        File_minmax <- tryCatch({
        
          # Load file
          raster <- rast(path)
        
          # Get min max
          File_minmax <- minmax(raster)
        
        }, error=function(e){File_minmax <- NA})
      
        cat("Finished calculating minmax for", path, "\n")
        return(File_minmax)
      }, simplify = FALSE)
    
      # Save minmax values for NCP
      saveRDS(NCP_minmax, NCP_minmax_path)
    }
  })
      
  # if report_NAs is TRUE, check for NAs and NaNs in the minmax values
  if(report_NAs == TRUE){
    
    # Summarise results across NCPs in a single dataframe
    All_NCP_minmaxs <- lapply(NCPs_to_summarise, function(NCP){
    
      # Load the minmax values for current NCP
      NCP_minmax <- readRDS(file.path(NCP_normalisation_dir, paste0(NCP, "_minmaxs.rds")))
  
      #loop over list and add details to dataframe
      NCP_minmax <- lapply(NCP_minmax, function(x){
        data.frame(Min = x[1], Max = x[2])
      })
  
      # bind list of dataframes into a single dataframe
      NCP_minmax <- as.data.frame(rbindlist(NCP_minmax, idcol = "Path"))
    })
    names(All_NCP_minmaxs) <- NCPs_to_summarise

    # bind list of dataframes into a single dataframe
    All_NCP_minmaxs <- as.data.frame(rbindlist(All_NCP_minmaxs, idcol = "NCP"))
    
    # Identify NAs
    NA_indices <- which(is.na(All_NCP_minmaxs$Min))

    # Get paths of the NA value
    NA_records <- NCP_layer_paths[NA_indices, "Path"]

    # Check for NaN values
    NaN_indices <- which((is.nan(All_NCP_minmaxs$Min)))

    # What is the path of the NaN value
    NaN_records <- NCP_layer_paths[NaN_indices, "Path"]
  return(c(NA_records, NaN_records))
  }
  }
  

#' calc_global_minmaxs:
#'
#' Calculate the global minmax values for each NCP
#'
#' This function calculates global minimum and maximum values
#' (i.e. smallest minimum and greatest maximum value) for each NCP across the
#' configurations*time points and saves the results to an RDS file
#'
#'
#' @param NCPs_to_summarise Vector of NCPs to summarise
#' @param NCP_normalisation_dir Directory to load/save the results
#'
#' @export

calc_global_minmaxs <- function(
    NCPs_to_summarise,
    NCP_normalisation_dir 
    ){
  
  # Outer loop over NCPs
  NCP_global_minmaxs <- lapply(NCPs_to_summarise, function(NCP){
    
    # Load the minmax values for current NCP
    NCP_minmax <- readRDS(file.path(NCP_normalisation_dir, paste0(NCP, "_minmaxs.rds")))
  
    #loop over list and add details to dataframe
    NCP_minmax <- lapply(NCP_minmax, function(x){
      data.frame(Min = x[1], Max = x[2])
    })
  
    # bind list of dataframes into a single dataframe
    NCP_minmax <- as.data.frame(rbindlist(NCP_minmax, idcol = "Path"))
    
    # Initialize min and max values
    NCP_min <- Inf
    NCP_max <- -Inf
  
    # Calculate global min and max for current NCP
    for(i in 1:nrow(NCP_minmax)){ 
        NCP_min <- min(NCP_min, NCP_minmax[i, "Min"], na.rm = TRUE)
        NCP_max <- max(NCP_max, NCP_minmax[i, "Max"], na.rm = TRUE)
    }
    return(data.frame(Min = NCP_min, Max = NCP_max))
  })
  names(NCP_global_minmaxs) <- NCPs_to_summarise

  # bind list of dataframes into a single dataframe
  NCP_global_minmaxs <- as.data.frame(rbindlist(NCP_global_minmaxs, idcol = "NCP"))

  # Save the global minmax values
  saveRDS(NCP_global_minmaxs, file = file.path(NCP_normalisation_dir, "NCP_global_minmaxs.rds"))
}

#' calc_layer_summaries:
#' 
#' Calculate the sum, mean and standard deviation of the normalised values
#'  for each time step of each configuration for the NCPs specificed 
#' 
#' This function calculates the sum, mean and standard deviation of the normalised values
#' for each time step of each configuration for the NCPs specificed and saves the results
#' to an RDS file
#' 
#' @param metrics Vector of metrics to calculate must use built-in function arguments from terra::global
#' i.e. "max", "min", "mean", "sum", "range", "rms" (root mean square), "sd",
#'  "std" (population sd, using n rather than n-1),
#' @param NCP_layer_paths Dataframe: containing info for each NCP layer
#' i.e. config_ID, time step, raw layer path and path to save normalized layer
#' @param NCPs_to_summarise Vector of NCPs to summarise
#' @param NCP_normalisation_dir Directory to load the global minmaxs values used for normalization 
#' @param NCP_summary_stats_dir Directory to save the summary statistics for individual NCPs
#' @param NCP_summarisation_dir Directory to save the summary statistics for all NCPs
#' @param Recalc_summary Logical indicating whether to recalculate the summary statistics for each NCP
#' @param Recalc_normalised_layers Logical indicating whether to recalculate the normalised layers
#' @param Save_normalised_layers Logical indicating whether to save the normalised layers
#' @param Sim_ctrl_tbl Dataframe: containing the simulation information
#' @param Normalise_results Logical indicating whether to normalise the results
#' @param Parallel Logical indicating whether to use parallel processing
#' 
#' @export
#

calc_layer_summaries <- function(
    metrics = c("sum", "mean", "sd"),
    NCP_layer_paths,
    NCPs_to_summarise,
    NCP_normalisation_dir,
    NCP_summary_stats_dir,
    NCP_summarisation_dir,
    Recalc_summary = TRUE,
    Recalc_normalised_layers = FALSE,
    Save_normalised_layers = FALSE,
    Sim_ctrl_tbl, 
    Normalise_results = TRUE,
    Parallel = TRUE
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
  NCP_global_minmaxs <- readRDS(file.path(NCP_normalisation_dir, "NCP_global_minmaxs.rds"))
  
  # Loop over the NCPs specified and calculate the summary metrics specified for the
  # normalized values for each time step under all configurations
  # for each configuration and time-step
  lapply(NCPs_to_summarise, function(NCP){
  
    # subset to NCP records
    NCP_layer_paths <- NCP_layer_paths[NCP_layer_paths$NCP == NCP,]
  
    # create path to save the combined results for this NCP
    NCP_sum_stats_path <- file.path(NCP_summary_stats_dir, paste0(NCP, "_summary_stats.rds"))
  
    # check if the summary file already exists & if Recalc_summary == FALSE then skip
    # otherwise if the file does not exist or Recalc_summary = TRUE then re-calculate
    if(file.exists(NCP_sum_stats_path) & Recalc_summary == FALSE){
      cat("Summary stats file already exists for", NCP, "\n")
      cat("Recalculation of summary is set to FALSE, skipping to next NCP \n")
    } else if (!file.exists(NCP_sum_stats_path) | Recalc_summary == TRUE){
    
      if(!file.exists(NCP_sum_stats_path)){
        cat("Summary stats file does not exist for", NCP, "\n")
      }
      if(Recalc_summary == TRUE){
        cat("Recalculation of summary is set to TRUE, calculating Summary stats for", NCP, "\n")
      }
    
      #if Recalc_normalised_layers is TRUE, then identify which layers exist
      #already and remove them from the df
      if(Recalc_normalised_layers == TRUE){

        cat("Recalculation of normalised layers is set to FALSE, identifying which layers already exist for", NCP, "\n")
    
        # Count how many of the normalized layers for this NCP exist and how many
        # are missing
        NCP_layer_paths$Exists <- sapply(NCP_layer_paths$Norm_path, file.exists)

        cat("Number of normalized raster layers for", NCP, "that exist:", sum(NCP_layer_paths$Exists), "\n")
        cat("Number of normalized raster layers for", NCP, "that are missing:", sum(!NCP_layer_paths$Exists), "\n") 
      
        # Remove the existing layers from the dataframe
        NCP_layer_paths <- NCP_layer_paths[!NCP_layer_paths$Exists,]
        
        cat("Existing layers have been removed from the list to be normalized, performing normalization of remaining layers \n")
      } else if(Recalc_normalised_layers == FALSE){
        cat("Recalculation of normalised layers is set to TRUE, normalising all layers in NCP_layer_paths \n")
      }
      
      # Loop over sequence of NCP layer paths, normalise the layer values 
      # with the global min/max values and then calculate sum, mean, median, mode
      # std of the normalised values. 
      NCP_sum_stats <- future_lapply(1:nrow(NCP_layer_paths), function(i){
  
        # use trycatch to handle missing or corrupt files
        Layer_info_stat <- tryCatch({
  
          # Load the raster
          raster <- rast(NCP_layer_paths[i, "Path"])
  
          NCP_name <- NCP_layer_paths[i, "NCP"]
          NCP_min <- NCP_global_minmaxs[NCP_global_minmaxs$NCP == NCP, "Min"]
          NCP_max <- NCP_global_minmaxs[NCP_global_minmaxs$NCP == NCP, "Max"]
  
          # Normalise the raster
          raster_norm <- (raster - NCP_min) / (NCP_max - NCP_min)
        
          # if Save_normalised_layers is TRUE, save the normalised raster
          if(Save_normalised_layers == TRUE){
            
            cat("Save_normalised_layers is set to TRUE \n")
            
            #create dir
            dir <- dirname(NCP_layer_paths[i, "Norm_path"])
            if(!dir.exists(dir)){
              dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            }

            # save as tif using the pre-prepared path
            writeRaster(raster_norm, file = NCP_layer_paths[i, "Norm_path"])
            
            cat("Normalized layer saved to", NCP_layer_paths[i, "Norm_path"], "\n")
          }
  
          # Calculate the stats
          raster_stats <- global(raster_norm, metrics, na.rm=TRUE)

          # cbind raster_stats with the row of NCP_layer_paths
          Layer_info_stat <- cbind(NCP_layer_paths[i,], raster_stats)
          }, error=function(e){

          # if an error occurs, return NA for all metrics  
          raster_stats <- rep(NA, length(metrics))
          names(raster_stats) <- metrics
          
          #convert to dataframe
          raster_stats <- t(data.frame(raster_stats))
          
          # cbind raster_stats with the row of NCP_layer_paths  
          Layer_info_stat <- cbind(NCP_layer_paths[i,], raster_stats)
        })
        return(Layer_info_stat)
      })
      
      # Save the results for the current NCP
      saveRDS(NCP_sum_stats, file = NCP_sum_stats_path)
      cat("Summary stats for", NCP, "have been saved to", NCP_sum_stats_path, "\n")
    }
  
  })

  # if parallel == TRUE set back to sequential processing
  if(Parallel == TRUE){
    future::plan(future::sequential)
  }

  # Produce a single df summarising the results across all NCPs
  NCP_sum_stats <- lapply(NCPs_to_summarise, function(NCP){
    
    # Load the summary stats for the current NCP
    NCP_sum_stats <- readRDS(file.path(NCP_summary_stats_dir, paste0(NCP, "_summary_stats.rds")))
    
    # Bind the list of dataframes for each time point together
    NCP_sum_stats <- do.call(rbind, NCP_sum_stats)
    
    return(NCP_sum_stats)
  })

  #rbind list of dataframes into a single dataframe
  NCP_sum_stats <- do.call(rbind, NCP_sum_stats)

  # Use the simulation control table to add the scenario_ID to the NCP_sum_stats
  NCP_sum_stats$Scenario <- sapply(NCP_sum_stats$Config_ID, function(x) Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == x, "Scenario_ID.string"])

  # Convert metric columns to numeric
  NCP_sum_stats[,metrics] <- lapply(NCP_sum_stats[,metrics], as.numeric)

  # Save the NCP_sum_stats
  saveRDS(NCP_sum_stats, file = file.path(NCP_summarisation_dir, "NCP_normalized_summary_stats.rds"))
  
  cat("Summary stats have been calculated and saved \n")
  
  # if Normalise_results is TRUE, then normalise the results and save a copy
  if(Normalise_results == TRUE){
    cat("Normalised_summary_stats is set to TRUE, normalising the summary stats \n")
    
    # Normalise the summary stats
    NCP_sum_stats_norm <- lapply(NCPs_to_summarise, function(NCP){
  
      # Get summary values for current NCP
      NCP_sum_stats_NCP <- NCP_sum_stats[NCP_sum_stats$NCP == NCP,]
  
      # Loop over the metrics
      for(metric in metrics){
    
         # Initialize min and max values
        NCP_min <- min(NCP_sum_stats_NCP[[metric]], na.rm = TRUE)
        NCP_max <- max(NCP_sum_stats_NCP[[metric]], na.rm = TRUE)
    
        # Normalize the summary values
        NCP_sum_stats_NCP[[paste0(metric, "_norm")]] <- (NCP_sum_stats_NCP[[metric]] - NCP_min) / (NCP_max - NCP_min)
  
        # Remove the original summary values
        NCP_sum_stats_NCP[[metric]] <- NULL
    
        # Rename the normalized summary values
        names(NCP_sum_stats_NCP)[names(NCP_sum_stats_NCP) == paste0(metric, "_norm")] <- metric
      }
  
    return(NCP_sum_stats_NCP)
    })
  names(NCP_sum_stats_norm) <- NCPs_to_summarise

  # bind the list of dataframes into a single dataframe
  NCP_sum_stats_norm <- as.data.frame(rbindlist(NCP_sum_stats_norm))

  # Save the normalised summary stats
  saveRDS(NCP_sum_stats_norm, file = file.path(NCP_summarisation_dir, "NCP_normalized_summary_stats_normalized.rds"))
  } 
}


#' calc_change_sumnmaries
#'
#' Calculate several metrics describing change in NCP provision between time
#' points:
#' 1. Avg. and std. of the positive and negative change in NCP 
#' provision for each NCP specified between each time point across all configurations
#' 2. a raster of the sum of change values across all time points each NCP and
#' all configurations
#' 3. The SSIM spatial pattern metric for simulated time point with reference
#' to the 1st time point for each NCP under all configurations: 
#' The SSIM calculation returns 4 metrics:
#' SIM: Similarity in Mean
#' SIV: Similarity in Variance
#' SIP: Similarity in Pattern of spatial covariance
#' SSIM: Which is the product of the three above
#' See research/Explanation_SSIM_components_Jones2016.png for details'
#'
#' @param NCP_layer_paths A dataframe with the paths to the NCP rasters for each
#' NCP and time point
#' @param NCPs_to_summarise A vector of NCPs to summarise
#' @param NCP_normalisation_dir The directory where the normalised NCP rasters are stored
#' @param NCP_change_summary_stats_dir The directory where the change summary stats are stored
#' @param NCP_change_sum_rasters_dir The directory where the change sum rasters are stored
#' @param NCP_SSIM_dir The directory where the SSIM results are stored
#' @param NCP_summarisation_dir The directory where the summarised results are stored
#' @param Sim_ctrl_tbl The simulation control table
#' @param Normalise_results A boolean indicating whether to normalise the summary stats
#' @param Parallel A boolean indicating whether to use parallel processing
#' @param Save_normalised_layers A boolean indicating whether normalised layers
#'  have been saved by the previous function calc_layer_summaries

calc_change_summaries <- function(
  NCP_layer_paths,
  NCPs_to_summarise,
  NCP_normalisation_dir,
  NCP_change_summary_stats_dir,
  NCP_change_sum_rasters_dir,
  NCP_SSIM_dir,
  NCP_summarisation_dir,
  Sim_ctrl_tbl, 
  Normalise_results = TRUE,
  Parallel = TRUE,
  Save_normalised_layers = TRUE
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
  
  # Load back in the global min maxs
  NCP_global_minmaxs <- readRDS(file.path(NCP_normalisation_dir, "NCP_global_minmaxs.rds"))
  
  # Outer loop over NCPs_to_summarise
  Calc_NCP_change <- lapply(NCPs_to_summarise, function(NCP){
  
    # Print message of current NCP
    cat(paste0("Calculating change summary stats for ", NCP, "\n"))
  
    # Get the global min and max for the NCP
    NCP_min <- NCP_global_minmaxs[NCP_global_minmaxs$NCP == NCP, "Min"]
    NCP_max <- NCP_global_minmaxs[NCP_global_minmaxs$NCP == NCP, "Max"]

    # Create the file paths for all expected outputs and check how many exist
    # Change summary stats    
    NCP_change_paths <- sapply(Config_IDs, function(Config_ID) file.exists(file.path(NCP_change_summary_stats_dir, paste0(NCP, "_Config_", Config_ID, "_change_stats.rds"))))
  
    # Print the number of files that exist
    cat(paste0("Number of configurations for which change summary stats have been computed for ", NCP, ": ", sum(NCP_change_paths), "\n"))
  
    # print number remaining
    cat(paste0("Number of configurations for which change summary stats remain to be computed for ", NCP, ": ", length(Config_IDs) - sum(NCP_change_paths), "\n"))
    
    # sum change raster
    NCP_sum_change_paths <- sapply(Config_IDs, function(Config_ID) file.exists(file.path(NCP_change_sum_rasters_dir, paste0(NCP, "_Config_", Config_ID, "_change_sum_rast.rds"))))
  
    # Print status of configs completed
    cat(paste0("Number of configurations for which sum change rasters have been computed for ", NCP, ": ", sum(NCP_sum_change_paths), "\n"))
    cat(paste0("Number of configurations for which change summary stats remain to be computed for ", NCP, ": ", length(Config_IDs) - sum(NCP_sum_change_paths), "\n"))
    
    # SSIM analysis
    NCP_SSIM_paths <- sapply(Config_IDs, function(Config_ID) file.exists(file.path(NCP_SSIM_dir, paste0(NCP, "_Config_", Config_ID, "_SSIM.rds"))))

    # Print status of configs completed
    cat(paste0("Number of configurations for which SSIM analysis has been computed for ", NCP, ": ", sum(NCP_SSIM_paths), "\n"))
    cat(paste0("Number of configurations for which SSIM analysis remains to be computed for ", NCP, ": ", length(Config_IDs) - sum(NCP_SSIM_paths), "\n"))
  
    # if all files exists then skip to next NCP
    if(all(NCP_change_paths) & all(NCP_sum_change_paths) & all(NCP_SSIM_paths)){
      cat(paste0("Change analysis and SSIM files exist for all configurations under ", NCP, ". Skipping to next NCP.\n"))
    } else {
    
      # Subset the Config_IDs for which change summary stats have not been computed
      # get indices of FALSE entries in NCP_change_paths, NCP_sum_change_paths and NCP_SSIM_paths
      Config_IDs_to_do <- unique(c(which(!NCP_change_paths), which(!NCP_sum_change_paths), which(!NCP_SSIM_paths)))
    
      # Sort in ascending order
      Config_IDs_to_do <- sort(Config_IDs_to_do)
    
      # Inner loop over Config_IDs
      Config_change_results <- future_lapply(Config_IDs_to_do, function(Config_ID){
    
        # Use tryCatch to catch errors and continue with the loop
        Config_results <- tryCatch({ 
      
        # Print message of current NCP and Config_ID  
        cat(paste0("Calculating change summary stats for ", NCP, " under Config ", Config_ID, "\n"))
    
        # Create the path for output of change analysis
        NCP_change_path <- file.path(NCP_change_summary_stats_dir, paste0(NCP, "_Config_", Config_ID, "_change_stats.rds"))
    
        # Create the path for the sum change raster
        NCP_sum_change_path <- file.path(NCP_change_sum_rasters_dir, paste0(NCP, "_Config_", Config_ID, "_change_sum_rast.rds"))
    
        # Create the path for output of SSIM analysis
        NCP_SSIM_path <- file.path(NCP_SSIM_dir, paste0(NCP, "_Config_", Config_ID, "_SSIM.rds"))
    
        # Check which of the files exist
        existing_files <- c(file.exists(NCP_change_path) & file.exists(NCP_sum_change_path) & file.exists(NCP_SSIM_path))

        # Check if all of the files exist if not then continue
        if(all(existing_files) == TRUE){
      
          # Print message that change analysis files already exist
          cat(paste0("All change and SSIM analysis files already exist for ", NCP, " under Config ", Config_ID, " skipping \n"))
        } else {
          
          # Get the paths for the NCP layers for the current NCP and Config_ID
          Config_info <- NCP_layer_paths[NCP_layer_paths$NCP == NCP & NCP_layer_paths$Config_ID == Config_ID,]
                    
          # if Save_normalised_layers is TRUE then normalized raster layers will
          # possibly have been saved already which means the normalisation 
          # does not need to be repeated (saving time)
          if(Save_normalised_layers == TRUE){
                    
            # Check if the normalised layers have been saved
            Config_info <- Config_info[unlist(lapply(Config_info$Norm_path, file.exists)),]  
                      
            # Read in the all the time step layers for the current NCP and Config_ID
            # as a list of rast objects
            NCP_layers <- lapply(Config_info$Norm_path, rast)
            
            if(nrow(Config_info) == 0){
              return(msg <- paste0("No input files found for Config_ID: ", Config_ID))
            }
          } else {
                      
            # Check that each of the time step files for this config exist and remove those that don't
            Config_info <- Config_info[unlist(lapply(Config_info$Path, file.exists)),]
                      
            if(nrow(Config_info) == 0){
            return(msg <- paste0("No input files found for Config_ID: ", Config_ID))
            }
            
            # Read in the all the time step layers for the current NCP and Config_ID
            # as a list of rast objects
            NCP_layers <- lapply(Config_info$Path, rast)
            
            # Normalize the layers
            NCP_layers <- lapply(NCP_layers, function(x){
              Layer_norm <- (x - NCP_min) / (NCP_max - NCP_min)
            })
          }
          
          # Change layer names to Time_steps
          names(NCP_layers) <- Config_info$Time_step
    
          # Check if both of the change analysis output files exist
          Change_files_exist <- c(file.exists(NCP_change_path), file.exists(NCP_sum_change_path))

          # If either file does not exist then calculate the change statistics
          if(all(Change_files_exist) == TRUE){
            cat(paste0("Change summary stats and sum raster already exist for ", NCP, " under Config ", Config_ID, " skipping \n"))
          } else {

            # Print message that change summary stats are being calculated
            cat(paste0("Either the change summary or the sum change raster does not exist, calculating change summary stats for ", NCP, " under Config ", Config_ID, "\n"))

            # Calculate the differences between each layer (i.e. differences over time)
            Diffs_over_time <- diff(rast(NCP_layers))

            # Calculate the average value of positive change cells
            Change_stats <- lapply(as.list(Diffs_over_time), function(Layer){
            vals <- values(Layer, na.rm = TRUE)
            avg_pos_change <- mean(vals[vals > 0])
            avg_neg_change <- mean(vals[vals < 0])
            std_pos_change <- sd(vals[vals > 0])
            std_neg_change <- sd(vals[vals < 0])

            # Return named list of results
            return(list(Avg_pos_change = avg_pos_change,
                        Avg_neg_change = avg_neg_change,
                        Std_pos_change = std_pos_change,
                        Std_neg_change = std_neg_change))
            })
            names(Change_stats) <- names(Diffs_over_time)

            # Bind to long dataframe
            Change_stats_long <- rbindlist(Change_stats, idcol = "Time_step")

            # Save as rds
            saveRDS(Change_stats_long, NCP_change_path)

            # Calculate sum of change over time points
            Change_sum <- sum(Diffs_over_time)

            # Save raster of the sum of differences over time as an RDS file
            saveRDS(Change_sum, NCP_sum_change_path)
            }
      
          #Check if the SSIM file exists and if not perform the calculation
          if(file.exists(NCP_SSIM_path) == TRUE){
            cat(paste0("SSIM file already exists for ", NCP, " under Config ", Config_ID, " skipping \n"))
          } else {

              # Print message that SSIM is being calculated
              cat(paste0("Calculating SSIM for ", NCP, " under Config ", Config_ID, "\n"))

              # Apply SSIM calculation to all layers using the 2020 (1st) layer as the reference
              SSIM_time_steps <- lapply(NCP_layers[2:length(NCP_layers)], function(x) {

                # Calculate SSIM
                SSIM_step <- capture.output(SSIMmap::ssim_raster(NCP_layers[[1]], x))

                # remove spaces in the output following ':'
                SSIM_step_str <- gsub(": ", ":", SSIM_step)

                # Split the output to a vector on the whitespace
                SSIM_step_str <- unlist(strsplit(SSIM_step_str, " "))

                # Name the list for the characters before the colon
                names(SSIM_step_str) <- sapply(SSIM_step_str, function(x) {
                  strsplit(x, ":")[[1]][1]
                })

                # split on the colon and only keep the numeric after it
                SSIM_step_str <- sapply(SSIM_step_str, function(x) {
                  as.numeric(strsplit(x, ":")[[1]][2])
                  })

                return(SSIM_step_str)
              })

              # bind results
              SSIM_time_steps_comb <- as.data.frame(do.call(rbind, SSIM_time_steps))

              #convert row names to a column
              SSIM_time_steps_comb$Time_step <- rownames(SSIM_time_steps_comb)

              # Save the result
              saveRDS(SSIM_time_steps_comb, NCP_SSIM_path)
            }

        

        # Print completion message
        cat(paste0("Completed change analysis for ", NCP, " under Config ", Config_ID, "\n"))
      }
    
   
      }, error=function(e){
        
        # Print error message
        cat(paste0("Error in Config ", Config_ID, ": ", e$message, "\n"))
    }) 
  })
  
    }
  })

  # if parallel == TRUE set back to sequential processing
  if(Parallel == TRUE){
    future::plan(future::sequential)
  }

  # Combine the change stats files together
  NCP_change_stats <- lapply(NCPs_to_summarise, function(NCP){
  
    # Create the file paths for all expected outputs and check how many exist
    NCP_change_paths <- lapply(Config_IDs, function(Config_ID) file.path(NCP_change_summary_stats_dir, paste0(NCP, "_Config_", Config_ID, "_change_stats.rds")))
  
    # Check if the change files exist
    NCP_change_paths_exist <- sapply(NCP_change_paths, function(x) file.exists(x))
    
    # Subset to only the existing files
    NCP_change_paths <- NCP_change_paths[NCP_change_paths_exist]
  
    # Subset Config IDs to only the existing files
    Config_IDs <- Config_IDs[NCP_change_paths_exist]
  
    # Load the existing files   
    NCP_change_dfs <- lapply(NCP_change_paths, function(x) readRDS(x))
    names(NCP_change_dfs) <- Config_IDs
  
    # bind results together
    NCP_change_summary <- rbindlist(NCP_change_dfs, idcol = "Config_ID")
    })
  names(NCP_change_stats) <- NCPs_to_summarise

  # bind results together
  NCP_change_stats <- as.data.frame(rbindlist(NCP_change_stats, idcol = "NCP"))

  # Add scenarios by matching on Config_ID
  NCP_change_stats$Scenario <- sapply(NCP_change_stats$Config_ID, function(x) Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == x, "Scenario_ID.string"])

  # Save the result
  saveRDS(NCP_change_stats, file.path(NCP_summarisation_dir, "NCPs_change_summary.rds"))

  # Combine the SSIM files together
  NCP_SSIM_stats <- lapply(NCPs_to_summarise, function(NCP){
  
    # Create the file paths for all expected outputs and check how many exist
    NCP_SSIM_paths <- lapply(Config_IDs, function(Config_ID) file.path(NCP_SSIM_dir, paste0(NCP, "_Config_", Config_ID, "_SSIM.rds")))
  
    # Check if the change files exist
    NCP_SSIM_paths_exist <- sapply(NCP_SSIM_paths, function(x) file.exists(x))
  
    # Subset to only the existing files
    NCP_SSIM_paths <- NCP_SSIM_paths[NCP_SSIM_paths_exist]
  
    # Subset Config IDs to only the existing files
    Config_IDs <- Config_IDs[NCP_SSIM_paths_exist]
  
    # Load the existing files   
    NCP_SSIM_dfs <- lapply(NCP_SSIM_paths, function(x) readRDS(x))
    names(NCP_SSIM_dfs) <- Config_IDs
  
    # bind results together
    NCP_SSIM_summary <- rbindlist(NCP_SSIM_dfs, idcol = "Config_ID")
    })
  names(NCP_SSIM_stats) <- NCPs_to_summarise

  # bind results together
  NCP_SSIM_stats <- as.data.frame(rbindlist(NCP_SSIM_stats, idcol = "NCP"))

  # Add scenarios by matching on Config_ID
  NCP_SSIM_stats$Scenario <- sapply(NCP_SSIM_stats$Config_ID, function(x) Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == x, "Scenario_ID.string"])
  
  # Save the result
  saveRDS(NCP_SSIM_stats, file.path(NCP_summarisation_dir, "NCPs_SSIM_summary.rds"))

  # if Normalise_results == TRUE, normalise the results
  if(Normalise_results == TRUE){
    
    # Get names of metrics from colnames removing NCP, Config_ID
    metrics_change <- colnames(NCP_change_stats)[!colnames(NCP_change_stats) %in% c("NCP", "Config_ID", "Scenario", "Time_step")]

    # Loop over NCPs_to_summarise and normalize the summary values of NCP provision
    NCP_change_stats_norm <- lapply(NCPs_to_summarise, function(NCP){
  
      # Get summary values for current NCP
      NCP_change_stats_NCP <- NCP_change_stats[NCP_change_stats$NCP == NCP,]
  
      # Loop over the metrics
      for(metric in metrics_change){
    
        # Initialize min and max values
        NCP_min <- min(NCP_change_stats_NCP[[metric]], na.rm = TRUE)
        NCP_max <- max(NCP_change_stats_NCP[[metric]], na.rm = TRUE)
    
        # for the Avg_neg_change metric, smaller values are better and hence we want to invert it so these have higher values
        if(metric == "Avg_neg_change"){
      
          # first change the values to absolute values
          NCP_change_stats_NCP[[metric]] <- abs(NCP_change_stats_NCP[[metric]])
      
          # Initialize min and max values
          NCP_min <- min(NCP_change_stats_NCP[[metric]], na.rm = TRUE)
          NCP_max <- max(NCP_change_stats_NCP[[metric]], na.rm = TRUE)
    
          # then invert the values
          NCP_change_stats_NCP[[metric]] <- (NCP_max - NCP_change_stats_NCP[[metric]] + NCP_min)
        }
    
        # Normalize the summary values
        NCP_change_stats_NCP[[paste0(metric, "_norm")]] <- (NCP_change_stats_NCP[[metric]] - NCP_min) / (NCP_max - NCP_min)
  
        # Remove the original summary values
        NCP_change_stats_NCP[[metric]] <- NULL
    
        # Rename the normalized summary values
        names(NCP_change_stats_NCP)[names(NCP_change_stats_NCP) == paste0(metric, "_norm")] <- metric
      }
  
      return(NCP_change_stats_NCP)
    })
    names(NCP_change_stats_norm) <- NCPs_to_summarise

    # bind the list of dataframes into a single dataframe
    NCP_change_stats_norm <- as.data.frame(rbindlist(NCP_change_stats_norm))

    # Save the NCP_change_stats_norm
    saveRDS(NCP_change_stats_norm, file = file.path(NCP_summarisation_dir, "NCP_normalized_change_stats_normalized.rds"))

    # Get names of metrics from colnames removing NCP, Config_ID
    metrics_SSIM <- colnames(NCP_SSIM_stats)[!colnames(NCP_SSIM_stats) %in% c("NCP", "Config_ID", "Scenario", "Time_step")]

    # Loop over NCPs_to_summarise and normalize the summary values of NCP provision
    NCP_SSIM_stats_norm <- lapply(NCPs_to_summarise, function(NCP){
  
      # Get summary values for current NCP
      NCP_SSIM_stats_NCP <- NCP_SSIM_stats[NCP_SSIM_stats$NCP == NCP,]
  
      # Loop over the metrics
      for(metric in metrics_SSIM){
    
        # Initialize min and max values
        NCP_min <- min(NCP_SSIM_stats_NCP[[metric]], na.rm = TRUE)
        NCP_max <- max(NCP_SSIM_stats_NCP[[metric]], na.rm = TRUE)
    
        # Normalize the summary values
        NCP_SSIM_stats_NCP[[paste0(metric, "_norm")]] <- (NCP_SSIM_stats_NCP[[metric]] - NCP_min) / (NCP_max - NCP_min)
  
        # Remove the original summary values
        NCP_SSIM_stats_NCP[[metric]] <- NULL
    
        # Rename the normalized summary values
        names(NCP_SSIM_stats_NCP)[names(NCP_SSIM_stats_NCP) == paste0(metric, "_norm")] <- metric
      }
  
  return(NCP_SSIM_stats_NCP)
})
    names(NCP_SSIM_stats_norm) <- NCPs_to_summarise

    # bind the list of dataframes into a single dataframe
    NCP_SSIM_stats_norm <- as.data.frame(rbindlist(NCP_SSIM_stats_norm))

    # Save the NCP_SSIM_stats_norm
    saveRDS(NCP_SSIM_stats_norm, file = file.path(NCP_summarisation_dir, "NCP_normalized_SSIM_normalized.rds"))
  }
}

#' Calc_aggregated_metric
#' 
#' The aggregated metric is calculated as the product of:
#' - Normalized SSIM to represent change in spatial pattern of NCP provision: 
#'     Rationale: Captures spatial dislocation of NCP provision from current state
#'     Interpretation: High values are good because they represent low spatial dislocation 
#' - Normalized average positive change in NCP provision: 
#'     Rationale: Capturing variation in disruption in provision
#'     Interpretation: High values are good because they represent higher average increases in provision
#' - Normalized average negative change in NCP provision
#'     Rationale: Capturing variation in disruption in provision
#'     Interpretation: High values are bad because they represent higher average decreases in provision NEED TO INVERT
#' - Normalized sum of NCP provision: 
#'     Rationale: Capturing total landscape level provision of NCPs
#'     Interpretation: High values are good because they represent higher total provision
#'     
#'  @param NCPs_to_summarise A character vector of NCPs to summarise
#'  @param NCP_summarisation_dir The directory where the NCP summarisation files are stored
#'  @param Agg_metrics A character vector of metrics to aggregate

calc_agg_metric <- function(
  NCPs_to_summarise,
  NCP_summarisation_dir,
  Agg_metrics = c("sum", "Avg_pos_change", "Avg_neg_change", "SSIM")
  ){
  
  # Load the NCP_sum_stats_norm
  Sum_stats_norm <- readRDS(file.path(NCP_summarisation_dir, "NCP_normalized_summary_stats_normalized.rds"))

  # Load the NCP_change_stats_norm
  Change_stats_norm <- readRDS(file.path(NCP_summarisation_dir, "NCP_normalized_change_stats_normalized.rds"))

  # Load the NCP_SSIM_stats_norm
  SSIM_stats_norm <- readRDS(file.path(NCP_summarisation_dir, "NCP_normalized_SSIM_normalized.rds"))

  # Loop over NCPs_to_summarise and calculate the aggregate metric
  NCP_Agg_stats <- lapply(NCPs_to_summarise, function(NCP){
  
    # Get summary values for current NCP
    Sum_stats_NCP <- Sum_stats_norm[Sum_stats_norm$NCP == NCP,]
    Change_stats_NCP <- Change_stats_norm[Change_stats_norm$NCP == NCP,]
    SSIM_stats_NCP <- SSIM_stats_norm[SSIM_stats_norm$NCP == NCP,]
  
    # create an empty list to store the seperate metrics
    NCP_Agg_stats_NCP <- list()
  
    # Loop over the metrics
    for(metric in Agg_metrics){
    
      # For the metrics in Sum_stats_norm
      if(metric %in% names(Sum_stats_norm)){
      
        # remove the other metrics from the dataframe and rename the metric column to be value
        NCP_stats <- Sum_stats_NCP %>% select(Config_ID, NCP, Scenario, Time_step, metric) %>% rename(value = metric)

        # Add the dataframe to the list
        NCP_Agg_stats_NCP[[metric]] <- NCP_stats
    
      # For the Change stats metrics
      } else if (metric %in% names(Change_stats_NCP)){ 

        # remove the other metrics from the dataframe
        NCP_stats <- Change_stats_NCP %>% select(Config_ID, NCP, Scenario, Time_step, metric)
      
        # if metric is Avg_neg_change then convert the negative values to positive
        if(metric == "Avg_neg_change"){
          NCP_stats$Avg_neg_change <- abs(NCP_stats[[metric]])
        }

        # rename the metric column to be value
        NCP_stats <- NCP_stats %>% rename(value = metric)
      
        # Add the dataframe to the list
        NCP_Agg_stats_NCP[[metric]] <- NCP_stats
      
      # For the SSIM metrics
      } else if (metric %in% names(SSIM_stats_NCP)){
      
        # remove the other metrics from the dataframe
        NCP_stats <- SSIM_stats_NCP %>% select(Config_ID, NCP, Time_step, Scenario, metric)
      
        # rename the metric column to be value
        NCP_stats <- NCP_stats %>% rename(value = metric)
      
        # Add the dataframe to the list
        NCP_Agg_stats_NCP[[metric]] <- NCP_stats
      }
      }
  
    #bind the lists together
    NCP_Agg_stats_NCP <- rbindlist(NCP_Agg_stats_NCP, idcol = "metric", use.names = TRUE)
  
    # pivot wide on the value columns with names from metric
    NCP_Agg_stats_NCP <- NCP_Agg_stats_NCP %>% pivot_wider(names_from = metric, values_from = value)
    
    # remove all entries of the 1st time point as the values will all be NA for this
    NCP_Agg_stats_NCP <- NCP_Agg_stats_NCP %>% filter(Time_step != min(Time_step))
    
    # loop over config_id and multiply all the metric columns together
    NCP_Agg_stats_NCP$Agg_metric <- sapply(1:nrow(NCP_Agg_stats_NCP), function(i){
      product <- prod(NCP_Agg_stats_NCP[i, Agg_metrics])
    })
  
    # get the min and max of the aggregated metric
    min_agg <- min(NCP_Agg_stats_NCP$Agg_metric, na.rm = TRUE)
    max_agg <- max(NCP_Agg_stats_NCP$Agg_metric, na.rm = TRUE)
  
    # Normalise the aggregated metric
    NCP_Agg_stats_NCP$Agg_metric_norm <- (NCP_Agg_stats_NCP$Agg_metric - min_agg)/ (max_agg - min_agg)
  
    return(NCP_Agg_stats_NCP)
  })
  names(NCP_Agg_stats) <- NCP_abbrevs

  # bind the list of dataframes into a single dataframe
  NCP_Agg_stats <- as.data.frame(rbindlist(NCP_Agg_stats))

  # Save the results
  saveRDS(NCP_Agg_stats, file.path(NCP_summarisation_dir, "NCP_aggregated_metric.rds"))
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
config <- config$Summarisation # only summarisation variables

# if InputDir is not set, use bash variable
# $FUTURE_EI_OUTPUT_DIR/$NCP_OUTPUT_BASE_DIR
if (is.null(config$InputDir) || config$InputDir == "") {
  config$InputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR,
                               bash_vars$NCP_OUTPUT_BASE_DIR)
}
# if OutputDir is not set, use bash variable SUMMARISATION_OUTPUT_DIR
if (is.null(config$OutputDir) || config$OutputDir == "") {
  config$OutputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR,
                                bash_vars$SUMMARISATION_OUTPUT_DIR)
}
# Check if InputDir is now set
if (is.null(config$InputDir) || config$InputDir == "") {
  stop("InputDir nor NCP_OUTPUT_BASE_DIR set in FUTURE_EI_CONFIG_FILE")
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

cat("Performing summarisation of Future EI outputs \n")

cat("Working directory is:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")

### Start by creating sub directories

# Dir for the processed NCP layers 
Normalized_NCP_layer_dir <- file.path(config$OutputDir, "Normalized_NCP_layers")
if (!dir.exists(Normalized_NCP_layer_dir)) {
  dir.create(Normalized_NCP_layer_dir, recursive = TRUE)
}

# Sub-dir for the summarisation results
NCP_summarisation_dir <- file.path(config$OutputDir, "NCP_summarisation")
if (!dir.exists(NCP_summarisation_dir)) {
  dir.create(NCP_summarisation_dir, recursive = TRUE)
}

# Sub-dir for the normalisation results
NCP_normalisation_dir <- file.path(NCP_summarisation_dir, "NCP_normalisation")
if (!dir.exists(NCP_normalisation_dir)) {
  dir.create(NCP_normalisation_dir, recursive = TRUE)
}

# Sub-dir for the summary_stats results
NCP_summary_stats_dir <- file.path(NCP_summarisation_dir, "NCP_summary_stats")
if (!dir.exists(NCP_summary_stats_dir)) {
  dir.create(NCP_summary_stats_dir, recursive = TRUE)
}

# Dir to store results of change analysis
NCP_change_dir <- file.path(NCP_summarisation_dir, "NCP_change")
if (!dir.exists(NCP_change_dir)) {
  dir.create(NCP_change_dir, recursive = TRUE)
}

# Sub-dir for Change_summary_stats and 
NCP_change_summary_stats_dir <- file.path(NCP_change_dir, "summary_stats")
if (!dir.exists(NCP_change_summary_stats_dir)) {
  dir.create(NCP_change_summary_stats_dir, recursive = TRUE)
}

# Sub-dir for rasters summing change over time
NCP_change_sum_rasters_dir <- file.path(NCP_change_dir, "sum_change_rasters")
if (!dir.exists(NCP_change_sum_rasters_dir)) {
  dir.create(NCP_change_sum_rasters_dir, recursive = TRUE)
}

# Sub-dir to store results of SSIM pattern analysis
NCP_SSIM_dir <- file.path(NCP_summarisation_dir, "NCP_SSIM")
if (!dir.exists(NCP_SSIM_dir)) {
  dir.create(NCP_SSIM_dir, recursive = TRUE)
}


### Prepare a dataframe of information on the NCP layers that are to be processed

# Load simulation control table
Sim_ctrl_tbl <- read.csv(config$Sim_control_table)

# Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI-")

# The Scenario_IDs from 865:1080 are actually BAU with climate change so rename 'BAU-CC'
Sim_ctrl_tbl$Scenario_ID.string[865:1080] <- "BAU-CC"

# Load EI interventions table
EI_interventions_tbl <- read.csv(config$EI_interventions_table)

# Get earliest scenario start date and latest end date
Start_date <- min(Sim_ctrl_tbl$Scenario_start.real)
End_date <- max(Sim_ctrl_tbl$Scenario_end.real)

# Create seq of scenario time steps with Step_length.real
Sim_time_steps <- seq(Start_date, End_date, by = Sim_ctrl_tbl$Step_length.real[1])

# Vector IDs of configurations to be analysed
# Note Manually adjust this to analyse specific configurations
Config_IDs <- unique(Sim_ctrl_tbl$Simulation_num.)

# Loop over NCPs_to_summarise and create vectors of file paths for their layers
# according to Config_IDs, Sim_time_steps and nesting
NCP_layer_paths <- lapply(config$NCPs_to_summarise, function(NCP){

  # Loop over config_IDs returning paths as vector
  Config_paths <- lapply(Config_IDs, function(Config_ID){

    # If NCP is not nested then append time step to file path
    if(config[["NCP_nesting"]][[NCP]] == FALSE){

      # Create a vector of file paths for each time step
      NCP_time_paths <- sapply(Sim_time_steps, function(Time_step){
        NCP_step_path <- file.path(config$InputDir, Config_ID, NCP, paste0(config[["NCP_file_names"]][[NCP]], "_", Time_step, ".tif"))
      })

    # If NCP is nesting then append time step as a dir before file path
    } else if (config[["NCP_nesting"]][[NCP]] == TRUE){

      # Create a vector of file paths for each time step
      NCP_time_paths <- sapply(Sim_time_steps, function(Time_step){

        NCP_step_path <-file.path(config$InputDir, Config_ID, NCP, Time_step, paste0(config[["NCP_file_names"]][[NCP]], ".tif"))

      })
    }

    # Also create a vector of file paths for saving the normalized NCP layers
    NCP_norm_paths <- sapply(Sim_time_steps, function(Time_step){
      file.path(Normalized_NCP_layer_dir , Config_ID, NCP, Time_step, paste0("Config_", Config_ID,"_", NCP, "_", Time_step, ".tif"))
    })
    
    # combine both vectors of paths in a dataframe
    NCP_time_paths <- data.frame(Path = NCP_time_paths, Norm_path = NCP_norm_paths)

    #Add column for NCP
    NCP_time_paths$NCP <- NCP

    # Add column for Config_ID
    NCP_time_paths$Config_ID <- Config_ID

    # Add column for Time_step
    NCP_time_paths$Time_step <- Sim_time_steps
  
    return(NCP_time_paths)
  })

  # Bind list of dataframes into a single dataframe
  NCP_df <- do.call(rbind, Config_paths)
})
names(NCP_layer_paths) <- config$NCPs_to_summarise

#rbind list of dataframes into a single dataframe
NCP_layer_paths <- do.call(rbind, NCP_layer_paths)

# Run calc_minmaxs function
calc_minmaxs(parallel = config$Parallel,
             NCP_layer_paths = NCP_layer_paths,
             NCP_normalisation_dir = NCP_normalisation_dir,
             NCPs_to_summarise = config$NCPs_to_summarise,
             minmax_recalc = TRUE, 
             report_NAs = TRUE)

# Run calc_global_minmaxs
calc_global_minmaxs(
    NCPs_to_summarise = config$NCPs_to_summarise,
    NCP_normalisation_dir = NCP_normalisation_dir)

# Run calc_layer_summaries
calc_layer_summaries(
    metrics = c("sum", "mean", "sd"),
    NCP_layer_paths = NCP_layer_paths,
    NCPs_to_summarise = config$NCPs_to_summarise,
    NCP_normalisation_dir = NCP_normalisation_dir,
    NCP_summary_stats_dir = NCP_summary_stats_dir,
    NCP_summarisation_dir = NCP_summarisation_dir,
    Recalc_summary = TRUE,
    Recalc_normalised_layers = FALSE,
    Save_normalised_layers = FALSE,
    Sim_ctrl_tbl = Sim_ctrl_tbl, 
    Normalise_results = TRUE,
    Parallel = config$Parallel
  )

# Run calc_change_summaries
calc_change_summaries(
  NCP_layer_paths = NCP_layer_paths,
  NCPs_to_summarise = config$NCPs_to_summarise,
  NCP_normalisation_dir = NCP_normalisation_dir,
  NCP_change_summary_stats_dir = NCP_change_summary_stats_dir,
  NCP_change_sum_rasters_dir = NCP_change_sum_rasters_dir,
  NCP_SSIM_dir = NCP_SSIM_dir,
  NCP_summarisation_dir = NCP_summarisation_dir,
  Sim_ctrl_tbl = Sim_ctrl_tbl, 
  Normalise_results = TRUE,
  Parallel = config$Parallel,
  Save_normalised_layers = FALSE
  )

# Run calc_agg_metric
calc_agg_metric(
  NCPs_to_summarise = config$NCPs_to_summarise,
  NCP_summarisation_dir = NCP_summarisation_dir,
  Agg_metrics = c("sum", "Avg_pos_change", "Avg_neg_change", "SSIM")
  )


