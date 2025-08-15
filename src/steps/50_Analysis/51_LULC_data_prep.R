
### =========================================================================
### Preparation
### =========================================================================

# Load libraries
packs <- c("stringr", "terra", "future", "future.apply", "readxl",
           "data.table", "tidyr", "yaml", "dplyr", "viridis", "ggplot2",
           "tidyterra", "treemapify", "jsonlite", "magick", "grDevices")
invisible(lapply(packs, require, character.only = TRUE))

# source save_indexed_png function
source("save_indexed_png.R")


web_platform_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform"
# if dir doesn't exist, create it
if(web_platform_dir){
  dir.create(web_platform_dir, recursive = TRUE, showWarnings = FALSE)
}

# Mask directory
Mask_dir <- file.path(web_platform_dir, "Masks")
# if dir doesn't exist, create it
if(!dir.exists(Mask_dir)){
  dir.create(Mask_dir, recursive = TRUE, showWarnings = FALSE)
}

### =========================================================================
### Prepare canton Bern mask
### =========================================================================

# Load shapefile of swiss cantons
canton_shp <- vect("Data/CH_geoms/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")

# filter to NAME == Bern
canton_bern <- canton_shp[canton_shp$NAME == "Bern", ]

# save the canton Bern shapefile to the mask directory
writeVector(canton_bern,
            file.path(Mask_dir, "Canton_mask.shp"),
             overwrite = TRUE)


### =========================================================================
### Finalising LULC tifs and images
### =========================================================================

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

prepare_lulc_files <- function(
    lulcc_input_dir = "F:/KB-outputs/lulcc_output",
    image_dir = "map_images",
    raster_dir = "raster_data",
    chart_data_dir = "chart_data",
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
  if(!dir.exists(file.path(base_dir, chart_data_dir))){
    dir.create(file.path(base_dir, chart_data_dir), recursive = TRUE, showWarnings = FALSE)
  }
  
  # load the simulation control table
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)
  
  # Tidy names in Sim_ctrl_tbl$Scenario_ID.string putting a '-' after 'EI
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "EI", "EI-")
  
  # for Scenario_ID.string replace any GREX with GR-EX
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(Sim_ctrl_tbl$Scenario_ID.string, "GREX", "GR-EX")
  
  # subset to only Scenario_ID.string == "GR-EX"
  Sim_ctrl_tbl <- Sim_ctrl_tbl[Sim_ctrl_tbl$Scenario_ID.string == "GR-EX", ]
  
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
  
  # Loop over config_IDs creating a dataframe with paths for saving the rasters, pngs and chart data
  lulc_paths <- lapply(Config_IDs, function(Config_ID){
    
    # Create a vector of file paths for each time step
    lulc_time_paths <- sapply(Sim_time_steps, function(Time_step){
      
      # Create the file path for the LULC layer: simulated_LULC_simID_179_year_2020
      lulc_step_path <- file.path(lulcc_input_dir, Config_ID, paste0( "simulated_LULC_simID_", Config_ID, "_year_", Time_step, ".tif"))
    })
    
    # use Config_ID to subset Sim_ctrl_tbl to get scenario details
    Config_details <- Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == Config_ID, ]
    
    # model path structure: NCP-Time_step-Config_details$Climate_scenario.string-Config_details$Econ_scenario.string-Config_details$Pop_scenario.string-Config_details$Scenario_ID.string-full.tif’
    lulc_paths_tif <- sapply(Sim_time_steps, function(Time_step){
      
      # ‘lulc-2020-rcp26-low-ref_central-bau-full.tif’
      file.path(base_dir, raster_dir,
                paste0("lulc-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                       Config_details$Econ_scenario.string, "-", 
                       Config_details$Pop_scenario.string, "-", 
                       Config_details$Scenario_ID.string, "-full.tif")) 
    })
    
    lulc_paths_png <- sapply(Sim_time_steps, function(Time_step){
      # ‘lulc-2020-rcp26-low-ref_central-bau-full.tif’
      file.path(base_dir, image_dir,
                paste0("lulc-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                       Config_details$Econ_scenario.string, "-", 
                       Config_details$Pop_scenario.string, "-", 
                       Config_details$Scenario_ID.string, "-full.png")) 
    })
    
    lulc_paths_perc_area <- sapply(Sim_time_steps, function(Time_step){
      # ‘lulc-2020-rcp26-low-ref_central-bau-full.tif’
      file.path(base_dir, chart_data_dir,
                paste0("lulc-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                       Config_details$Econ_scenario.string, "-", 
                       Config_details$Pop_scenario.string, "-", 
                       Config_details$Scenario_ID.string, "-full-perc_area.json")) 
    })
    
    lulc_paths_area_chg <- sapply(Sim_time_steps, function(Time_step){
      # ‘lulc-2020-rcp26-low-ref_central-bau-full.tif’
      file.path(base_dir, chart_data_dir,
                paste0("lulc-", Time_step, "-", Config_details$Climate_scenario.string, "-", 
                       Config_details$Econ_scenario.string, "-", 
                       Config_details$Pop_scenario.string, "-", 
                       Config_details$Scenario_ID.string, "-full-area_chg.json")) 
      
    })
    
    # combine all vectors of paths in a dataframe
    lulc_time_paths <- data.frame(Path = lulc_time_paths,
                                  lulc_path_tif = lulc_paths_tif, 
                                  lulc_path_png = lulc_paths_png,
                                  lulc_path_perc_area = lulc_paths_perc_area,
                                  lulc_path_area_chg = lulc_paths_area_chg
    )
    
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
    
    #load scenario specific glacier index
    Glacier_index <- readRDS(file = list.files("Data/glacier_scenario_indices",
                                               full.names = TRUE,
                                               pattern = Scenario))[,c("ID_loc", paste(Start_date))]
    
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
              
              # modify the save paths for the masked raster layer, png and json files
              mask_path_tif <- gsub("full", mask_name, scenario_lulc_df$lulc_path_tif[i])
              map_path_png <- gsub("full", mask_name, scenario_lulc_df$lulc_path_png[i])
              mask_path_data <- gsub("full", mask_name, scenario_lulc_df$lulc_path_perc_area[i])
              
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
              
              # apply the mask to the raster layer
              masked_lulc <- terra::mask(x = lulc_layer, 
                                         mask = mask_layer, 
                                         updatevalue = NA)
              
              
              # update the file path for the masked raster layer by replacing 'full'
              # in scenario_lulc_df$lulc_path_tif[i] with the mask name
             
              
              # save the masked raster layer
              writeRaster(masked_lulc,
                          filename = mask_path_tif,
                          overwrite = TRUE)
              
              cat(paste("Saved masked raster layer to:", mask_path_tif, "\n"))
              
              # get frequency table of the masked raster layer
              rast_tbl <- freq(masked_lulc)
              rast_tbl$layer <- NULL
              rast_tbl$class_name <- c(unlist(sapply(rast_tbl$value, function(y) unique(unlist(Aggregation_scheme[Aggregation_scheme$Aggregated_ID == y, "Aggregated_class_short"])),simplify = TRUE)), "Lake", "River")
              rast_tbl$value <- NULL
              rast_tbl$perc_area <- rast_tbl$count / sum(rast_tbl$count) * 100
              rast_tbl$count <- NULL
              
              # modify the data path
              mask_path_data <- gsub("full", mask_name, scenario_lulc_df$lulc_path_perc_area[i])
              
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
    
    
    
  })
}

# test function
prepare_lulc_files(
  lulcc_input_dir = "F:/KB-outputs/lulcc_output",
  image_dir = "map_images",
  raster_dir = "raster_data",
  chart_data_dir = "chart_data",
  base_dir = web_platform_dir,
  Sim_ctrl_tbl_path = "LULCC_CH_HPC/simulation_control.csv",
  ProjCH = ProjCH,
  LULC_agg_path = "LULCC_CH_HPC/Tools/LULC_class_aggregation.xlsx",
  colour_pal = LULC_pal,
  Non_agg_lulc_path = "Data/NOAS04_2018.tif",
  Use_parallel = FALSE,
  num_workers = 4,
  map_masks = list("canton" = file.path(Mask_dir, "Canton_mask.shp")),
  overwrite = FALSE
)















    


