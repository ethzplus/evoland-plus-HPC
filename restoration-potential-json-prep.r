# Load libraries
packs <- c("stringr", "terra", "future", "future.apply", "readxl",
           "data.table", "tidyr", "yaml", "dplyr", "viridis", "ggplot2",
           "tidyterra", "jsonlite", "magick", "grDevices", "classInt", "doparallel", "foreach")
invisible(lapply(packs, require, character.only = TRUE))

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
                             save_samples,
                             break_types = c("quantile", "fisher", "regular"),
                             save_plots,
                             parallel,
                             input_map_path_name,
                             output_tbl_path_name,
                             regular_min = -1,
                             regular_max = 1,
                             regular_by = 0.25
                             ){
  
  # break_types includes quantiles or fisher then it is necessary to take a 
  # sample from each raster layer in order to calculate global breaks under these methods
  if(any(c("quantile", "fisher") %in% break_types)){
    if(parallel == TRUE){
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl)
      
      samples_list <- foreach(i = 1:nrow(ES_layer_paths), .combine = c, .packages = "terra") %dopar% {
        if(save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])){
          samp <- readRDS(ES_layer_paths[i, "sample_path"])
          return(samp)
        } else if(file.exists(ES_layer_paths[i, "sample_path"]) == FALSE){
          raster_norm <- rast(ES_layer_paths[i, "tif_path"])
          vals <- values(raster_norm, mat=FALSE)
          vals <- vals[!is.na(vals)]
          samp <- sample(vals, min(sample_size, length(vals)), replace=FALSE)
          if(save_samples == TRUE){
            saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
          }
          return(samp)
        }
      }
      stopCluster(cl)
      all_samples <- samples_list
    } else{
      all_samples <- c()
      for(i in 1:nrow(ES_layer_paths)){
        if(save_samples == TRUE & file.exists(ES_layer_paths[i, "sample_path"])){
          samp <- readRDS(ES_layer_paths[i, "sample_path"])
          all_samples <- c(all_samples, samp)
          next
        } else if(file.exists(ES_layer_paths[i, "sample_path"]) == FALSE){
          raster_norm <- rast(ES_layer_paths[i, input_map_path_name])
          vals <- values(raster_norm, mat=FALSE)
          vals <- vals[!is.na(vals)]
          samp <- sample(vals, min(sample_size, length(vals)), replace=FALSE)
          if(save_samples == TRUE){
            saveRDS(samp, file = ES_layer_paths[i, "sample_path"])
          }
          all_samples <- c(all_samples, samp)
        }
      }
    }
  }
  
  # default deciles for quantile/fisher
  probs <- seq(0, 1, 0.1) 
  
  # empty breaks list
  breaks <- list()
  for(break_type in break_types){
    breaks[[break_type]] <- NA
  }

  for(break_type in break_types){
    if(break_type == "quantile"){
      breaks$quantile <- unique(quantile(all_samples, probs=probs, na.rm=TRUE))
    } else if(break_type == "fisher"){
      fisher <- classIntervals(all_samples, n = length(probs), style = "fisher")
      breaks$fisher <- fisher$brks 
    } else if(break_type == "regular"){
      breaks$regular <- seq(regular_min, regular_max, by = regular_by)
    }
  }
  
  future_lapply(1:nrow(ES_layer_paths), function(i){
    raster <- rast(ES_layer_paths[i, input_map_path_name])
    
    for(break_type in names(breaks)){
      current_breaks <- breaks[[break_type]]
      
      if(break_type %in% c("quantile", "fisher")){
        rcl_matrix <- cbind(current_breaks[-length(current_breaks)], 
                            current_breaks[-1], 
                            1:(length(current_breaks)-1))
        r_classified <- terra::classify(raster, rcl = rcl_matrix)
      } else if(break_type == "regular"){
        rcl_matrix <- cbind(current_breaks[-length(current_breaks)], 
                            current_breaks[-1], 
                            1:(length(current_breaks)-1))
        r_classified <- terra::classify(raster, rcl = rcl_matrix)
      }
      
      freq_tbl <- freq(r_classified)
      freq_tbl <- as.data.frame(freq_tbl)
      freq_tbl$value[is.na(freq_tbl$value)] <- 0
      
      # remove value of 0
      freq_tbl <- freq_tbl[freq_tbl$value != 0, ]
      
      # if freq_tbl does not contain all unique values from column 3 in rcl_matrix
      # then add rows for the missing values with count 0 and layer ==1
      missing_values <- setdiff(1:(length(current_breaks)-1), freq_tbl$value)
      if(length(missing_values) > 0){
        missing_df <- data.frame(value = missing_values, count = 0, layer = 1)
        freq_tbl <- rbind(freq_tbl, missing_df)
      }
      
      
      if(nrow(freq_tbl) == 0) {
        cat("Warning: No data for", break_type, "- skipping\n")
        next
      }
      total_cells <- sum(freq_tbl$count)
      freq_tbl$percent <- 100 * freq_tbl$count / total_cells
      
      n_expected_bins <- length(current_breaks) - 1
      all_labels <- character(n_expected_bins)
      for(j in 1:n_expected_bins) {
        lower <- current_breaks[j]
        upper <- current_breaks[j + 1]
        all_labels[j] <- paste0(round(lower, 2), " – ", round(upper, 2))
      }
      
      freq_tbl$bin_label <- sapply(freq_tbl$value, function(bin_num) {
        if(is.na(bin_num)) {
          return("NA")
        } else if(bin_num < 1 || bin_num > length(all_labels)) {
          return(paste0("Bin_", bin_num, "_OutOfRange"))
        } else {
          return(all_labels[bin_num])
        }
      })

      freq_tbl$method <- break_type
      json_data <- toJSON(setNames(as.list(freq_tbl$bin_label), freq_tbl$percent), pretty = TRUE)
      json_path <- gsub(".json", paste0("_", break_type, ".json"), 
                        ES_layer_paths[i, output_tbl_path_name])
      write(json_data, file = ES_layer_paths[i, output_tbl_path_name])
      
      # sum of percent should be 100
      sum <- sum(freq_tbl$percent)
      if(abs(sum - 100) > 0.01){
        cat("Warning: Percentages do not sum to 100 for", ES_layer_paths[i, "ES"], "using", break_type, "breaks\n")
      }
      
      if(save_plots == TRUE){
        chart_images_dir <- gsub("chart_data", "chart_images", dirname(ES_layer_paths[i, output_tbl_path_name]))
        if(!dir.exists(chart_images_dir)){
          dir.create(chart_images_dir, recursive = TRUE, showWarnings = FALSE)
        }
        chart_path <- file.path(chart_images_dir, gsub(".json", ".png", basename(json_path)))
        bar_plot <- ggplot(freq_tbl, aes(x = bin_label, y = percent)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          labs(title = paste("Area in Classes for", ES_layer_paths[i, "ES"], "using", break_type, "breaks"),
               x = "Class",
               y = "Percentage of Area") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(filename = chart_path, plot = bar_plot, width = 10, height = 6)
      }
    }
  })
  cat("Finished calculating breaks for all layers of current ES and saving classified area JSON files.\n")
}

#' save_continuous_indexed_png
#' Save a continuous indexed PNG from a raster object with specified colors
#' called by normalise_layers function
#' @param raster_obj A SpatRaster object to be saved as PNG
#' @param output_path The file path where the PNG will be saved
#' @param low_color The color for the lowest value in the raster
#' @param high_color The color for the highest value in the raster
#' @param mid_color An optional color for the mid-point in the raster
#' @param mask An optional SpatRaster, SpatVector, file path to a raster, or file path to a shapefile defining the spatial mask. NA values within this mask will appear white, while NA values outside will be transparent
#' @param na_color The color to use for NA values within the mask (default is "white")
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
                                       mask = NULL,
                                       na_color = "white",
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
                                       verbose = FALSE) {
  
  # Load required libraries
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required but not installed.")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.")
  }
  
  library(magick)
  library(terra)
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
  
  # Process raster with mask if provided
  if (!is.null(mask)) {
    if (verbose) cat("Processing mask input...\n")
    
    # Load mask from file path if it's a character string
    if (is.character(mask)) {
      if (!file.exists(mask)) {
        stop(paste("Mask file does not exist:", mask))
      }
      
      # Determine file type and load accordingly
      file_ext <- tolower(tools::file_ext(mask))
      
      if (file_ext %in% c("tif", "tiff", "grd", "nc", "img")) {
        # Load as raster
        if (verbose) cat("Loading mask as raster...\n")
        mask <- rast(mask)
        
      } else if (file_ext == "shp") {
        # Load shapefile and rasterize
        if (verbose) cat("Loading shapefile and rasterizing to match input raster...\n")
        
        # Read shapefile
        shp <- vect(mask)
        
        # Transform CRS if needed
        if (crs(shp) != crs(raster_obj)) {
          if (verbose) cat("Transforming shapefile CRS to match raster...\n")
          shp <- project(shp, crs(raster_obj))
        }
        
        # Rasterize the shapefile to match the input raster
        mask <- rasterize(shp, raster_obj, field = 1)
        
      } else {
        stop(paste("Unsupported mask file format. Supported formats: .tif, .tiff, .grd, .nc, .img (rasters) and .shp (shapefiles). Got:", file_ext))
      }
    }
    
    # Ensure mask is a SpatRaster
    if (!inherits(mask, "SpatRaster")) {
      stop("Mask must be a SpatRaster, file path to a raster, or file path to a shapefile")
    }
    
    # Use only the first layer if multi-layer
    if (nlyr(mask) > 1) {
      if (verbose) cat("Using first layer of multi-layer mask...\n")
      mask <- mask[[1]]
    }
    
    # Ensure mask and raster have same extent and resolution
    if (!compareGeom(raster_obj, mask, lyrs = FALSE, crs = TRUE, warncrs = FALSE, ext = TRUE, rowcol = TRUE, res = TRUE)) {
      if (verbose) cat("Resampling mask to match raster...\n")
      mask <- resample(mask, raster_obj)
    }
    
    # Create a copy of the raster for processing
    processed_raster <- raster_obj
    
    # Find NA values within the mask (where mask has non-NA values)
    na_within_mask <- is.na(raster_obj) & !is.na(mask)
    
    # Assign a special value to NA pixels within the mask
    # We'll use a value outside the normal range to distinguish it
    special_na_value <- high_value + (high_value - low_value) * 0.1
    processed_raster[na_within_mask] <- special_na_value
    
    # Build color palette including the special NA color
    if (!is.null(mid_color)) {
      if (verbose) cat("Generating continuous palette with mid color and mask NA color...\n")
      # Reserve one color for the special NA value
      data_colors <- colorRampPalette(c(low_color, mid_color, high_color))(max_colors - 1)
      color_palette <- c(data_colors, na_color)
    } else {
      if (verbose) cat("Generating continuous palette with mask NA color...\n")
      # Reserve one color for the special NA value
      data_colors <- colorRampPalette(c(low_color, high_color))(max_colors - 1)
      color_palette <- c(data_colors, na_color)
    }
    
    # Update the zlim to include the special NA value
    plot_zlim <- c(low_value, special_na_value)
    plot_raster <- processed_raster
    
  } else {
    # No mask provided, use original behavior
    if (!is.null(mid_color)) {
      if (verbose) cat("Generating continuous palette with mid color...\n")
      color_palette <- colorRampPalette(c(low_color, mid_color, high_color))(max_colors)
    } else {
      if (verbose) cat("Generating continuous palette without mid color...\n")
      color_palette <- colorRampPalette(c(low_color, high_color))(max_colors)
    }
    
    plot_zlim <- c(low_value, high_value)
    plot_raster <- raster_obj
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
  
  # Plot with explicit zlim to ensure consistent color mapping
  plot(plot_raster, 
       col = color_palette, 
       range = plot_zlim,  # terra uses 'range' instead of 'zlim'
       legend = show_legend, 
       axes = axes, 
       box = box)
  
  dev.off()
  
  # Quantize to indexed PNG
  if (verbose) cat("Converting to indexed PNG...\n")
  img <- image_read(temp_png)
  img_indexed <- img %>%
    image_quantize(max = max_colors, colorspace = colorspace, dither = FALSE) %>%
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





# base dir for web platform
web_platform_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform"

# restoration potential dir
RP_dir <- file.path(web_platform_dir, "restoration_potential")

# list of restoration potential tifs
RP_file_paths <- list.files(file.path(RP_dir, "restoration_potential_tifs"), pattern = ".tif", full.names = TRUE)

# replace '.tif' with .json and use RP_dir for the base dir
RP_output_tbl_paths <- gsub(".tif", ".json", RP_file_paths)
RP_output_tbl_paths <- sapply(RP_output_tbl_paths, function(x) {
  file.path(RP_dir, basename(x))
})

RP_png_path <- gsub(".tif", ".png", RP_file_paths)
RP_png_path <- sapply(RP_png_path, function(x){
  file.path(RP_dir, basename(x))
})

# create a df called ES_layer_paths with columns: input_path and output_tbl_path
ES_layer_paths <- data.frame(
  tif_path = RP_file_paths,
  output_tbl_path = RP_output_tbl_paths,
  png_path = RP_png_path,
  stringsAsFactors = FALSE
)

# apply function
Calculate_areas_in_classes(ES_layer_paths,
                             sample_size = 50000,
                             save_samples = FALSE,
                             break_types = c("regular"),
                             save_plots = FALSE,
                             parallel = FALSE,
                             input_map_path_name = "tif_path",
                             output_tbl_path_name = "output_tbl_path",
                             regular_min = 0,
                             regular_max = 1,
                             regular_by = 0.1
                             )


# Mask directory
Mask_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/Masks"

# apply function to save images
for(i in 1:nrow(ES_layer_paths)){
  raster_obj <- rast(ES_layer_paths[i, "tif_path"])
  output_path <- ES_layer_paths[i, "png_path"]

  save_continuous_indexed_png(raster_obj, 
                              output_path, 
                              low_color = "#F2C7A6",
                              high_color = "#00A9A4",
                              mid_color = NULL,
                              low_value = 0,      # New argument for low value
                              high_value = 1,     # New argument for high value
                              mid_value = NULL,    # New argument for mid value (used when mid_color is specified)
                              mask = file.path(Mask_dir, "Canton_mask.shp"),        
                              na_color = "white",
                              width = 25, 
                              height = 20, 
                              resolution = 300,
                              units = "cm",
                              margins = c(0, 0, 0, 0),
                              background = "transparent",
                              #colorspace = "sRGB",
                              max_colors = 256,
                              show_legend = FALSE,
                              axes = FALSE,
                              box = FALSE,
                              cleanup_temp = TRUE,
                              verbose = FALSE)
   cat(paste("Saved PNG for", ES_layer_paths[i, "tif_path"], "to", output_path, "\n"))
}

