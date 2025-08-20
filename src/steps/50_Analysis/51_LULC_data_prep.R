#' Produce LULC maps and images for the web platform
#'
#' Input: Folder of outputs from LULCC simulations
#' Output: rasters, images and chart data for the web platform
#'
#' @environment ncp_summarisation
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
Sys.setenv(LULCC_CH_HPC_DIR = "LULCC_CH_HPC")
Sys.setenv(LULCC_M_EI_INTS_TABLE = "EI_interventions.csv")

# Load libraries
packs <- c("stringr", "terra", "future", "future.apply", "readxl",
           "data.table", "tidyr", "yaml", "dplyr", "viridis", "ggplot2",
           "tidyterra", "treemapify", "jsonlite", "magick", "grDevices")
invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

### =========================================================================
### Functions
### =========================================================================

#' Save Raster as 8-bit Indexed PNG
#'
#' This function saves a raster object as an 8-bit indexed color PNG file
#' using the magick package while preserving the original colors.
#'
#' @param raster_obj A raster object to be saved
#' @param output_path Character string specifying the output PNG file path
#' @param color_palette Vector of colors to use for the raster (hex colors)
#' @param width Numeric, width of the output image in specified units (default: 25)
#' @param height Numeric, height of the output image in specified units (default: 20)
#' @param resolution Numeric, resolution in DPI (default: 300)
#' @param units Character, units for width and height ("cm", "in", "px") (default: "cm")
#' @param margins Numeric vector of length 4, margins in the order c(bottom, left, top, right) (default: c(0, 0, 0, 0))
#' @param background Character, background color (default: "transparent")
#' @param colorspace Character, colorspace for quantization ("sRGB" or "rgb") (default: "sRGB")
#' @param max_colors Numeric, maximum number of colors for indexed palette (default: 256)
#' @param show_legend Logical, whether to show the raster legend (default: FALSE)
#' @param axes Logical, whether to show axes (default: FALSE)
#' @param box Logical, whether to show a box around the plot (default: FALSE)
#' @param cleanup_temp Logical, whether to remove temporary files (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#'
#' @return Invisibly returns the magick image object of the final indexed PNG
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' save_indexed_png(my_raster, "output.png", my_colors)
#' 
#' # With custom dimensions and margins
#' save_indexed_png(my_raster, "output.png", my_colors, 
#'                  width = 30, height = 25, margins = c(1, 1, 1, 1))
#' 
#' # High resolution output
#' save_indexed_png(my_raster, "output.png", my_colors, 
#'                  resolution = 600, units = "in", width = 8, height = 6)
#' }
#'
#' @export
save_indexed_png <- function(raster_obj, 
                           output_path, 
                           color_palette,
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
  
  if (verbose) cat("Creating temporary PNG...\n")
  
  # Create temporary PNG file
  temp_png <- tempfile(fileext = ".png")
  
  # Create the standard PNG first
  png(temp_png, 
      width = width, 
      height = height, 
      res = resolution, 
      bg = background, 
      units = units)
  
  # Set margins
  par(mar = margins)
  
  # Plot the raster
  plot(raster_obj, 
       col = color_palette, 
       legend = show_legend, 
       axes = axes, 
       box = box)
  
  dev.off()
  
  if (verbose) cat("Converting to indexed PNG...\n")
  
  # Read the temporary PNG with magick
  img <- image_read(temp_png)
  
  # Convert to 8-bit indexed color PNG while preserving colors
  img_indexed <- img %>%
    image_quantize(max = max_colors, 
                   colorspace = colorspace, 
                   dither = FALSE) %>%
    image_strip()  # Remove unnecessary metadata
  
  # Write the indexed PNG
  image_write(img_indexed, output_path, format = "png", depth = 8)
  
  if (verbose) cat(paste("Saved indexed PNG to:", output_path, "\n"))
  
  # Clean up temporary file
  if (cleanup_temp) {
    unlink(temp_png)
    if (verbose) cat("Cleaned up temporary files.\n")
  } else {
    if (verbose) cat(paste("Temporary file saved at:", temp_png, "\n"))
  }
  
  # Return the magick image object invisibly
  invisible(img_indexed)
}

# Helper function to inspect the saved PNG properties
#' Inspect Indexed PNG Properties
#' This function inspects the properties of an indexed PNG file
#' using both the magick and png packages.
#' @param png_path Character string specifying the path to the PNG file
#' @return Invisibly returns a list containing information from both magick and png packages
inspect_indexed_png <- function(png_path) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required but not installed.")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' is required but not installed.")
  }
  
  library(magick)
  library(png)
  
  cat("=== PNG Inspection Results ===\n")
  
  # Method 1: Using magick
  img <- image_read(png_path)
  info <- image_info(img)
  
  cat("Magick Information:\n")
  cat(paste("  Dimensions:", info$width, "x", info$height, "\n"))
  cat(paste("  Color space:", info$colorspace, "\n"))
  cat(paste("  Depth:", info$depth, "bit\n"))
  cat(paste("  Format:", info$format, "\n"))
  
  # Method 2: Using png package for detailed info
  png_info <- readPNG(png_path, info = TRUE)
  info_attr <- attr(png_info, "info")
  
  color_types <- c("0" = "Grayscale", 
                   "2" = "RGB", 
                   "3" = "Indexed", 
                   "4" = "Grayscale + Alpha", 
                   "6" = "RGBA")
  
  cat("\nDetailed PNG Information:\n")
  cat(paste("  Color Type:", color_types[as.character(info_attr$color.type)], 
            "(", info_attr$color.type, ")\n"))
  cat(paste("  Is Indexed:", info_attr$color.type == 3, "\n"))
  cat(paste("  Bit Depth:", info_attr$bit.depth, "\n"))
  
  return(invisible(list(magick_info = info, png_info = info_attr)))
}


#' Prepare LULC files for the web platform
#' This function prepares LULC files for the web platform by processing simulation outputs, creating rasters, images, and chart data.
#' @param lulcc_input_dir Directory containing LULCC simulation outputs
#' @param image_dir Directory to save map images
#' @param raster_dir Directory to save raster data
#' @param perc_area_data_dir Directory to save percentage area data
#' @param area_chg_data_dir Directory to save area change data
#' @param base_dir Base directory for the web platform
#' @param Sim_ctrl_tbl_path Path to the simulation control table CSV file
#' @param ProjCH Projection string for the raster data description
#' @param LULC_agg_path Path to the LULC aggregation scheme Excel file
#' @param colour_pal Named list of colors for LULC classes
#' @param Non_agg_lulc_path Path to the non-aggregated LULC raster file
#' @param Use_parallel Logical, whether to use parallel processing (default: FALSE)
#' @param num_workers Number of workers for parallel processing (default: 4)
#' @param map_masks List of masks to apply to the LULC maps (default: full Canton mask)
#' @param overwrite Logical, whether to overwrite existing files (default: TRUE)
#' 
prepare_lulc_files <- function(
    lulcc_input_dir = "F:/KB-outputs/lulcc_output",
    image_dir = "map_images",
    raster_dir = "raster_data",
    perc_area_data_dir = "chart_data/perc_area",
    area_chg_data_dir = "chart_data/perc_area_change",
    base_dir = web_platform_dir,
    Sim_ctrl_tbl_path = Sim_ctrl_tbl_path,
    ProjCH = ProjCH,
    LULC_agg_path = "LULCC_CH_HPC/Tools/LULC_class_aggregation.xlsx",
    colour_pal = LULC_pal,
    Non_agg_lulc_path = "Data/NOAS04_2018.tif",
    Use_parallel = FALSE,
    num_workers = 4,
    map_masks = list("full" = file.path(Mask_dir, "Canton_mask.shp")),
    overwrite = TRUE
){
  
  # if the directories do not exist, create them
  if(!dir.exists(file.path(base_dir, image_dir))){
    dir.create(file.path(base_dir, image_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, raster_dir))){
    dir.create(file.path(base_dir, raster_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, perc_area_data_dir))){
    dir.create(file.path(base_dir, perc_area_data_dir), recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(file.path(base_dir, area_chg_data_dir))){
    dir.create(file.path(base_dir, area_chg_data_dir), recursive = TRUE, showWarnings = FALSE)
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
  
  # Load the LULC aggregation scheme
  Aggregation_scheme <- read_excel(LULC_agg_path)
  
  # convert the values in Class_abbreviation to lower case
  Aggregation_scheme$Class_abbreviation <- tolower(Aggregation_scheme$Class_abbreviation)
  
  # Loop over config_IDs creating a dataframe with paths for saving the rasters, pngs and chart data
  lulc_paths <- lapply(Config_IDs, function(Config_ID){
    
    # use Config_ID to subset Sim_ctrl_tbl to get scenario details
      Config_details <- Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == Config_ID, ]
    
    # Create a dataframe of file paths for each time step
    lulc_time_paths <- as.data.frame(rbindlist(lapply(Sim_time_steps, function(Time_step){
      
      # Create the file path for the LULC layer: simulated_LULC_simID_179_year_2020
      input_path <- file.path(lulcc_input_dir, Config_ID, paste0( "simulated_LULC_simID_", Config_ID, "_year_", Time_step, ".tif"))
    
      # model path structure: NCP-Time_step-Config_details$Climate_scenario.string-Config_details$Econ_scenario.string-Config_details$Pop_scenario.string-Config_details$Scenario_ID.string-full.tif’
      tif_path <- file.path(base_dir, raster_dir,
                            paste0("lulc-", Time_step, "-",
                                   Config_details$Climate_scenario.string, "-", 
                                   Config_details$Econ_scenario.string, "-", 
                                   Config_details$Pop_scenario.string, "-", 
                                    Config_details$Scenario_ID.string, "-full.tif"))
      
      # ‘lulc-2020-rcp26-low-ref_central-bau-full.tif’
      png_path <- file.path(base_dir, image_dir,
                            paste0("lulc-", Time_step, "-",
                                  Config_details$Climate_scenario.string, "-", 
                                  Config_details$Econ_scenario.string, "-", 
                                  Config_details$Pop_scenario.string, "-", 
                                  Config_details$Scenario_ID.string, "-full.png"))
      
      perc_area_path <- file.path(base_dir, perc_area_data_dir,
                            paste0("lulc-", Time_step, "-", 
                                   Config_details$Climate_scenario.string, "-", 
                                   Config_details$Econ_scenario.string, "-", 
                                   Config_details$Pop_scenario.string, "-", 
                                   Config_details$Scenario_ID.string, "-full-perc_area.json"))
      
      # loop over the mask names in map_masks modifying the tif, png and perc_area paths
      mask_paths_tif <- sapply(names(map_masks), function(mask_name){
        gsub("full", mask_name, tif_path)
      })
      
      mask_paths_png <- sapply(names(map_masks), function(mask_name){
        gsub("full", mask_name, png_path)
      })
          
      mask_paths_perc_area <- sapply(names(map_masks), function(mask_name){
        gsub("full", mask_name, perc_area_path)
      })
        
      # combine all paths into a named list
      time_step_paths <- list(input_path, tif_path, png_path, perc_area_path, 
                                  mask_paths_tif, mask_paths_png, mask_paths_perc_area)
      names(time_step_paths) <- c("Path", "lulc_path_tif", "lulc_path_png", "lulc_path_perc_area",
                             paste0("lulc_path_tif_", names(map_masks)),
                             paste0("lulc_path_png_", names(map_masks)),
                             paste0("lulc_path_perc_area_", names(map_masks)))
      
      return(time_step_paths)
      })))

    # Add column for Config_ID
    lulc_time_paths$Config_ID <- Config_ID
    
    # Add column for Time_step
    lulc_time_paths$Time_step <- Sim_time_steps
    
    # add column for Scenario_ID
    lulc_time_paths$Scenario_ID <- Config_details$Scenario_ID.string
    
    # add column for Climate_scenario
    lulc_time_paths$Climate_scenario <- Config_details$Climate_scenario.string
    
    # add column for Econ_scenario
    lulc_time_paths$Econ_scenario <- Config_details$Econ_scenario.string
    
    # add column for Pop_scenario
    lulc_time_paths$Pop_scenario <- Config_details$Pop_scenario.string
    
    # loop over NCP_time_paths and check which exist
    lulc_time_paths$Exists <- file.exists(lulc_time_paths$Path)
    
    return(lulc_time_paths)
  })
  
  # Bind list of dataframes into a single dataframe
  lulc_df <- do.call(rbind, lulc_paths)
  
  # check if all files exist
  if (all(lulc_df$Exists)) {
    message("All LULC files exist.")
  } else {
    stop("Some LULC files do not exist. Please check the paths.")
  }
  
  # next step is to replace the glacier , rivers and lakes values in the LULC maps
  # Using one of the LULC rasters and the aggregation scheme to create a raster attribute table
  # get unique class values from raster and add 20 and 21 to represent Lake and River
  LULC_rat <- data.frame(ID = c(sort(unique(values(rast(lulc_df$Path[1])))), 20, 21))
  LULC_rat$lulc_name <- c(unlist(sapply(LULC_rat$ID, function(y) unique(unlist(Aggregation_scheme[Aggregation_scheme$Aggregated_ID == y, "Aggregated_class_short"])),simplify = TRUE)), "Lake", "River")  
  
  #subset aggregation table to distinct value of the aggregated LULC classes
  subset_agg <- Aggregation_scheme %>% distinct(Aggregated_ID, .keep_all=TRUE)
  
  #add colours to class info
  subset_agg$colours <- sapply(subset_agg$Aggregated_class_short, function(x){
    colour <- colour_pal[[paste(x)]]
  })
  
  LULC_rat$colour <- sapply(LULC_rat$lulc_name, function(x){
    colour <- colour_pal[[paste(x)]]
  })
  
  # Create a named vector for color mapping
  col_map <- setNames(LULC_rat$colour, LULC_rat$ID)
  
  #Load in most recent non-aggregated LULC raster
  ref_LULC <- rast(Non_agg_lulc_path)
  
  # check if the raster has the correct projection
  if (crs(ref_LULC) != ProjCH) {
    ref_LULC <- project(ref_LULC, ProjCH)
  }
  
  # get raster values of lakes and rivers
  mask_values <- unlist(Aggregation_scheme[Aggregation_scheme$NOAS04_class_ENG %in% c("Lakes", "Rivers", "Glaciers, perpetual snow"), "NOAS04_ID"])
  names(mask_values) <- c(20,21, 19)
  
  #convert raster to dataframe
  ref_LULC_dat <- terra::as.data.frame(ref_LULC, na.rm = FALSE) 
  
  #add ID column to dataset
  ref_LULC_dat$ID <- seq.int(nrow(ref_LULC))
  
  #Get XY coordinates of cells
  xy_coordinates <- crds(ref_LULC, na.rm=FALSE) 
  
  #cbind XY coordinates to dataframe and seperate rows where all values = NA
  ref_LULC_dat <- cbind(ref_LULC_dat, xy_coordinates)
  
  #vector the pixel values of Glacier from the non_agg lulc layer
  Glacier_val <- unlist(Aggregation_scheme[Aggregation_scheme$NOAS04_class_ENG == "Glaciers, perpetual snow", "NOAS04_ID"])
  
  # because we need to adjust the glacier locations in the reference lulc map
  #according to each scenario for efficiency use an outer loop over the scenarios
  # and an inner loop over the scenario specific simulations
  
  if(Use_parallel){
    plan(multisession, workers = num_workers)
  } else {
    plan(sequential)
  }
  
  future_sapply(unique(lulc_df$Scenario_ID), function(Scenario){
    
    cat(paste("Processing scenario:", Scenario, "\n"))
    
    Glacier_files <- list.files("Data/glacier_scenario_indices",
                                               full.names = TRUE)
    
    # match on the scenario ignoring the case
    Scenario_glacier_file <- Glacier_files[grepl(Scenario, Glacier_files, ignore.case = TRUE)]
  
    #load scenario specific glacier index
    Glacier_index <- readRDS(Scenario_glacier_file)[,c("ID_loc", paste(Start_date))]
    
    #seperate vector of cell IDs for glacier and non-glacer cells
    Non_glacier_IDs <- Glacier_index[Glacier_index[[paste(Start_date)]]==0, "ID_loc"]
    Glacier_IDs <- Glacier_index[Glacier_index[[paste(Start_date)]]==1, "ID_loc"] 
    
    # create a copy of ref_LULC_dat for this scenario
    scenario_LULC_dat <- ref_LULC_dat
    
    #replace the 1's and 0's with the correct LULC
    scenario_LULC_dat[scenario_LULC_dat$ID %in% Non_glacier_IDs, "NOAS04_2018"] <- unlist(Aggregation_scheme[Aggregation_scheme$NOAS04_class_ENG == "Rocks", "NOAS04_ID"])
    scenario_LULC_dat[scenario_LULC_dat$ID %in% Glacier_IDs, "NOAS04_2018"] <- Glacier_val
    
    #2nd step ensure that other glacial cells that do not match the glacier index
    #are also changed to static so that the transition rates calculate the
    #correct number of cell changes
    scenario_LULC_dat[which(scenario_LULC_dat$NOAS04_2018 == Glacier_val & !(scenario_LULC_dat$ID %in% Glacier_IDs)), "NOAS04_2018"] <- unlist(Aggregation_scheme[Aggregation_scheme$NOAS04_class_ENG == "Rocks", "NOAS04_ID"])   
    
    #convert back to raster
    scenario_LULC <- rast(scenario_LULC_dat[,c("x", "y", "NOAS04_2018")], crs = ProjCH)
    
    # seperate the scenario specific lulc_df
    scenario_lulc_df <- lulc_df[lulc_df$Scenario_ID == Scenario, ]
    
    # loop over the rows of lulc_df for this scenario
    for(i in 1:nrow(scenario_lulc_df)){
      
      # print the simulation id and time step
      cat(paste("Processing simulation:", scenario_lulc_df$Config_ID[i], 
                "at time step:", scenario_lulc_df$Time_step[i], "\n"))
      
      # check if the tif file exists
      if (scenario_lulc_df$Exists[i]) {
        
        # read the raster layer
        lulc_layer <- rast(scenario_lulc_df$Path[i])
        
        # add the crs to the raster layer
        crs(lulc_layer) <- ProjCH
        
        # if the time step is the first one, we need to replace the glacier, river and lake values
        # which are the current mask values if not then we only need to replace rivers and lakes
        if(scenario_lulc_df$Time_step[i] == min(scenario_lulc_df$Time_step)){
          time_step_mask_values <- mask_values
        } else {
          time_step_mask_values <- mask_values[names(mask_values) %in% c(20, 21)]
        }
        
        #loop over lulc mask values
        for(j in 1:length(mask_values)){
          lulc_layer <- terra::mask(x = lulc_layer,
                                    mask = scenario_LULC,
                                    maskvalues = time_step_mask_values[j],
                                    updatevalue = as.numeric(names(time_step_mask_values)[j]))
          
        } #close for loop over mask values
        
        # Now loop over any masks provided in the map_masks list
        if(length(map_masks) > 0){
          
          cat("Masking map to specificed areas \n")
          
          for(mask_name in names(map_masks)){
            
            cat(paste("Applying mask:", mask_name, "\n"))
            
            # Check if the mask file exists
            mask_path <- map_masks[[mask_name]]
            
            # If the mask file exists, apply it to the raster layer
            if(file.exists(mask_path)){
              
              # get the mask paths
              mask_path_tif <- scenario_lulc_df[i,paste0("lulc_path_tif_", mask_name)]
              map_path_png <- scenario_lulc_df[i,paste0("lulc_path_png_", mask_name)]
              mask_path_data <- scenario_lulc_df[i, paste0("lulc_path_perc_area_", mask_name)]
              
              # if all of these files already exist and overwrite == FALSE, skip to the next iteration
              if(!overwrite && 
                 file.exists(mask_path_tif) && 
                 file.exists(map_path_png) && 
                 file.exists(mask_path_data)){
                message(paste("Files already exist, skipping:", mask_path_tif, map_path_png, mask_path_data))
                next
              }
              
              # if the mask path contains shp extension, read it as a vector
              if(grepl("\\.shp$", mask_path)){
                message(paste("Applying mask from shapefile:", mask_path))
                mask_layer <- vect(mask_path)
              } else if(grepl("\\.tif$", mask_path)){
                message(paste("Applying mask from raster file:", mask_path))
                mask_layer <- rast(mask_path)
              } else {
                stop(paste("Unsupported mask file type for:", mask_path))
              }
              
              # crop the lulc_layer to the extent of the mask layer
              lulc_layer <- terra::crop(x = lulc_layer, 
                                        y = mask_layer)
              
              # apply the mask to the raster layer
              masked_lulc <- terra::mask(x = lulc_layer, 
                                         mask = mask_layer, 
                                         updatevalue = NA)
              
              # save the masked raster layer
              writeRaster(masked_lulc,
                          filename = mask_path_tif,
                          overwrite = TRUE)
              
              cat(paste("Saved masked raster layer to:", mask_path_tif, "\n"))
              
              # get frequency table of the masked raster layer
              rast_tbl <- freq(masked_lulc)
              rast_tbl$layer <- NULL
              rast_tbl$class_name <- c(unlist(sapply(rast_tbl$value, function(y) unique(unlist(Aggregation_scheme[Aggregation_scheme$Aggregated_ID == y, "Class_abbreviation"])),simplify = TRUE)), "lake", "river")
              rast_tbl$value <- NULL
              rast_tbl$perc_area <- rast_tbl$count / sum(rast_tbl$count) * 100
              rast_tbl$count <- NULL
              
              # convert the df to json
              json_data <- toJSON(setNames(as.list(rast_tbl$class_name), rast_tbl$perc_area), pretty = TRUE)
              
              # save the json data to the mask_path_data
              write(json_data, file = mask_path_data)
              
              cat(paste("Saved table of LULC % areas to:", mask_path_data, "\n"))
              
              # plot the raster layer matching the values to colours in LULC_rat
              
              # use the save_indexed_png function to save the raster layer as a png
              save_indexed_png(
                raster_obj = masked_lulc, 
                output_path = map_path_png, 
                color_palette = col_map,
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
              
              cat(paste("Saved map image to:", map_path_png, "\n"))
  
            } else {
              message(paste("Mask file does not exist:", mask_path))
            }
          }
        } else {
          message("No masks provided.")
        }
      } else {
        message(paste("File does not exist:", scenario_lulc_df$Path[i]))
      }
    }
    
    # now that we have process all time points for each individual simulation
    # we can calculate changes in the % of class areas
    # now that we have calculated the % area of each class for each time step
    # we can calculate the difference in % class area between simulation start and end
  
    # loop over each unique configuration ID for this scenario
    for(Config_ID in unique(scenario_lulc_df$Config_ID)){
      
      # subset the data to this Config_ID
      config_lulc_df <- scenario_lulc_df[scenario_lulc_df$Config_ID == Config_ID, ]
      
      cat(paste("Calculating area change for Config_ID:", Config_ID, "\n"))
      
      # get the first time step for this scenario
      first_time_step <-  min(config_lulc_df$Time_step)
      
      # list time steps for this Config_ID
      config_time_steps <- unique(config_lulc_df$Time_step)
      
      # remove the first time step from the list of time steps
      config_time_steps <- config_time_steps[config_time_steps != first_time_step]
      
      # loop over the time steps calculating the perc area change compared to the first time step for each
      for(time_step in config_time_steps){
        
        cat(paste("Calculating area change from", first_time_step, "to", time_step, "\n"))
        
        # loop over mask names in the map_masks list
        for(mask_name in names(map_masks)){
          
          # create the area change path for this mask name
          area_change_path <- file.path(base_dir, area_chg_data_dir, 
                                        paste0("lulc-",time_step, "-",
                                               config_lulc_df$Climate_scenario[config_lulc_df$Time_step == time_step], 
                                               "-",config_lulc_df$Econ_scenario[config_lulc_df$Time_step == time_step], 
                                               "-", config_lulc_df$Pop_scenario[config_lulc_df$Time_step == time_step], 
                                               "-", config_lulc_df$Scenario_ID[config_lulc_df$Time_step == time_step], 
                                               "-", mask_name, 
                                               "-perc_area_chg.json"))
          
          # if this path already exists and overwrite is FALSE, skip to the next iteration
          if(!overwrite && file.exists(area_change_path)){
            message(paste("Area change file already exists, skipping:", area_change_path))
            next
          }
          
          # check that there is a corresponding column for this mask name in the config_lulc_df
          if(!paste0("lulc_path_perc_area_", mask_name) %in% colnames(config_lulc_df)){
            stop(paste("No column for lulc_path_perc_area_", mask_name, "in config_lulc_df"))
          }
          
          # get the file path for the first step and current step
          first_mask_path <- config_lulc_df[config_lulc_df$Time_step == first_time_step, paste0("lulc_path_perc_area_", mask_name)]
          current_mask_path <- config_lulc_df[config_lulc_df$Time_step == time_step, paste0("lulc_path_perc_area_", mask_name)]
          
          # read the first and current time step data
          first_data <- fromJSON(first_mask_path, simplifyDataFrame = TRUE)
          
          # convert to df with names as a column called perc_area and values as a column called class_name
          first_data <- data.frame(class_name = unlist(first_data), perc_area = as.numeric(names(first_data)), stringsAsFactors = FALSE, row.names = NULL)
          
          current_data <- fromJSON(current_mask_path)
          current_data <- data.frame(class_name = unlist(current_data), perc_area = as.numeric(names(current_data)), stringsAsFactors = FALSE, row.names = NULL)
          
          
          # calculate the difference in % area for each class as a % of the area at the start
          area_change <- merge(first_data, current_data, by = "class_name", suffixes = c("_start", "_end"))
          
          
          area_change$perc_area_change <- ((area_change$perc_area_end - area_change$perc_area_start)/area_change$perc_area_start)*100
          
          # remove the lakes and rivers from the area change data
          area_change <- area_change[!area_change$class_name %in% c("Lake", "River"), ]
          
          # sort from largest to smallest change
          area_change <- area_change %>% arrange(desc(perc_area_change))
          
          # make sure the class_name is a factor with levels in the current order
          area_change$class_name <- factor(area_change$class_name, levels = area_change$class_name)
        
          # convert the df to json
          area_change_json <- toJSON(setNames(as.list(area_change$class_name), area_change$perc_area_change), pretty = TRUE)
              
          # save the json data to the mask_path_data
          write(area_change_json, file = area_change_path)
          
          cat(paste("Saved area change data to:", area_change_path, "\n"))
        } # close for loop over mask names
        } # close for loop over time steps
      
        cat("finished all time steps for Config_ID:", Config_ID, "\n")
    } # close for loop over unique Config_IDs
  })
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


# dir for saving results
web_platform_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform"
# if dir doesn't exist, create it
if(!dir.exists(web_platform_dir)){
  dir.create(web_platform_dir, recursive = TRUE, showWarnings = FALSE)
}

# Mask directory
Mask_dir <- file.path(web_platform_dir, "Masks")
# if dir doesn't exist, create it
if(!dir.exists(Mask_dir)){
  dir.create(Mask_dir, recursive = TRUE, showWarnings = FALSE)
}

#crs for maps
ProjCH <- "+proj=somerc +init=epsg:2056"

#define colour palette as list using LULC_rat$lulc_name as names
LULC_pal <- list("Urban/amenities" = '#a8aba5', #Urban
            "Static" = "#d1d3cf", #static
            "Open Forest" = "#97d1d5", #Open forest
            "Closed forest" = "#29898f", #closed forest
            "Overgrown/shrubland" = "#bb8a75", #Shrubland
            "Intensive agriculture" =  "#f59f78", #Intensive agriculture
            "Alpine pastures" = "#6ca147", #Alpine pastures
            "Grassland or meadows" = "#c4e0a1", #Grassland
            "Permanent crops" = "#DDCC66", #Permanet crops
            "Glacier" = "#d5f1ff",
            "River" = "#93d0ee",
            "Lake" = "#93d0ee")

### =========================================================================
### Prepare canton Bern mask
### =========================================================================

# Load shapefile of swiss cantons
canton_shp <- vect("Data/CH_geoms/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")

# filter to NAME == Bern
canton_bern <- canton_shp[canton_shp$NAME == "Bern", ]

# check to see that the CRS matches ProjCH
if (crs(canton_bern) != ProjCH) {
  # if not, transform the CRS
  canton_bern <- terra::project(canton_bern, ProjCH)
}

# save the canton Bern shapefile to the mask directory
writeVector(canton_bern,
            file.path(Mask_dir, "Canton_mask.shp"),
             overwrite = TRUE)


### =========================================================================
### Finalising LULC tifs and images
### =========================================================================

# Apply function to prepare LULC files
prepare_lulc_files(
    lulcc_input_dir = "F:/KB-outputs/lulcc_output",
    image_dir = "map_images",
    raster_dir = "raster_data",
    perc_area_data_dir = "chart_data/perc_area",
    area_chg_data_dir = "chart_data/perc_area_change",
    base_dir = web_platform_dir,
    Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
    ProjCH = ProjCH,
    LULC_agg_path = "LULC_class_aggregation.xlsx",
    colour_pal = LULC_pal,
    Non_agg_lulc_path = "Data/NOAS04_2018.tif",
    Use_parallel = FALSE,
    num_workers = 4,
    map_masks = list("canton" = file.path(Mask_dir, "Canton_mask.shp")),
    overwrite = TRUE
)















    


