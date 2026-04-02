#' Produce LULC maps and images for the web platform
#'
#' Input: Folder of outputs from LULCC simulations
#' Output: rasters, images and chart data for the web platform
#'
#' @environment ncp_summarisation
#' @config $FUTURE_EI_CONFIG_FILE (yaml file)
#' @date 2026-03-27
#' @author Benjamin Black
#'

# FOR TESTING ONLY
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
  "treemapify",
  "jsonlite",
  "magick",
  "grDevices"
)
invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

### =========================================================================
### Functions
### =========================================================================

#' Ensure a directory exists.
#' @param dir Directory path to create.
#' @return Invisibly returns `dir`.
ensure_dir <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  invisible(dir)
}

#' Save a raster as a PNG using a fixed value-to-colour mapping.
#' @param raster_obj A SpatRaster object to be saved as a PNG.
#' @param output_path The file path where the PNG should be saved.
#' @param color_palette Named vector mapping raster values to colours.
#' @param width The width of the output PNG in the specified units (default: 25).
#' @param height The height of the output PNG in the specified units (default: 20).
#' @param resolution The resolution of the output PNG in dots per inch (default: 300).
#' @param units The units for width and height (default: "cm").
#' @param margins A numeric vector of length 4 specifying the margins (bottom, left, top, right) in lines (default: c(0, 0, 0, 0)).
#' @param background The background colour for the PNG (default: "transparent").
#' @param colorspace The colour space for the PNG (default: "sRGB").
#' @param max_colors The maximum number of colours to use in the PNG (default: 256).
#' @param show_legend Logical indicating whether to include a legend in the plot (default: FALSE).
#' @param axes Logical indicating whether to include axes in the plot (default: FALSE).
#' @param box Logical indicating whether to include a box around the plot (default: FALSE).
#' @param cleanup_temp Logical indicating whether to delete the temporary PNG.
#' @param verbose Logical indicating whether to print verbose messages during processing (default: FALSE).
#' @return An invisible magick image object representing the saved PNG.
save_indexed_png <- function(
  raster_obj,
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
  verbose = FALSE
) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required but not installed.")
  }

  raster_values <- sort(unique(terra::values(raster_obj)))
  raster_values <- raster_values[!is.na(raster_values)]

  if (length(raster_values) == 0) {
    stop("Raster contains no non-NA values to plot.")
  }

  palette_keys <- as.character(raster_values)
  missing_keys <- setdiff(palette_keys, names(color_palette))
  if (length(missing_keys) > 0) {
    stop(
      paste(
        "Missing colour mapping for raster values:",
        paste(missing_keys, collapse = ", ")
      )
    )
  }

  ordered_palette <- unname(color_palette[palette_keys])

  if (verbose) {
    cat("Creating PNG with palette ordered by raster values\n")
  }

  temp_png <- tempfile(fileext = ".png")

  png(
    temp_png,
    width = width,
    height = height,
    res = resolution,
    bg = background,
    units = units
  )

  par(mar = margins)
  plot(
    raster_obj,
    col = ordered_palette,
    legend = show_legend,
    axes = axes,
    box = box
  )

  dev.off()

  img <- magick::image_read(temp_png)

  # Preserve the plotted colours exactly rather than re-quantizing and
  # risking palette reordering inside ImageMagick.
  img <- magick::image_strip(img)
  magick::image_write(img, output_path, format = "png", depth = 8)

  if (cleanup_temp) {
    unlink(temp_png)
  }

  invisible(img)
}

#' Compute percentage area by LULC class for a raster.
#' @param raster A SpatRaster object representing the LULC data.
#' @param Aggregation_scheme Data frame with `Aggregated_ID` and `Class_abbreviation`.
#' @param context_label Optional label used in log messages.
#' @return Named list of class names to percentage area, or `NULL` if no data.
process_lulc_area <- function(raster, Aggregation_scheme, context_label = "") {
  cat("Computing frequency table for", context_label, "\n")

  rast_tbl <- freq(raster)
  rast_tbl$layer <- NULL

  if (nrow(rast_tbl) == 0) {
    cat("Warning: No data in raster for", context_label, "- skipping\n")
    return(NULL)
  }

  class_names <- character(nrow(rast_tbl))

  for (i in seq_len(nrow(rast_tbl))) {
    val <- rast_tbl$value[i]

    if (val == 20) {
      class_names[i] <- "Lake"
    } else if (val == 21) {
      class_names[i] <- "River"
    } else {
      matched_class <- unique(unlist(Aggregation_scheme[
        Aggregation_scheme$Aggregated_ID == val,
        "Class_abbreviation"
      ]))

      if (length(matched_class) > 0) {
        class_names[i] <- matched_class[1]
      } else {
        class_names[i] <- paste0("Unknown_", val)
        cat("Warning: No class found for value", val, "in", context_label, "\n")
      }
    }
  }

  rast_tbl$class_name <- class_names
  rast_tbl <- rast_tbl[order(rast_tbl$value), ]
  rast_tbl$class_name <- factor(
    rast_tbl$class_name,
    levels = rast_tbl$class_name
  )

  rast_tbl$value <- NULL
  rast_tbl$perc_area <- rast_tbl$count / sum(rast_tbl$count) * 100
  rast_tbl$count <- NULL

  setNames(as.list(rast_tbl$perc_area), rast_tbl$class_name)
}

#' Load and normalise simulation control and aggregation inputs.
#' @param Sim_ctrl_tbl_path The file path to the simulation control table CSV file.
#' @param LULC_agg_path The file path to the LULC aggregation scheme Excel file.
#' @return List with control table, aggregation scheme, config IDs, time steps,
#'   and overall start/end dates.
load_and_prepare_config <- function(Sim_ctrl_tbl_path, LULC_agg_path) {
  Sim_ctrl_tbl <- read.csv(Sim_ctrl_tbl_path, stringsAsFactors = FALSE)

  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "EI",
    "EI_"
  )
  Sim_ctrl_tbl$Scenario_ID.string <- str_replace_all(
    Sim_ctrl_tbl$Scenario_ID.string,
    "GREX",
    "GR_EX"
  )

  Sim_ctrl_tbl$Climate_scenario.string <- tolower(
    Sim_ctrl_tbl$Climate_scenario.string
  )
  Sim_ctrl_tbl$Econ_scenario.string <- tolower(
    Sim_ctrl_tbl$Econ_scenario.string
  )
  Sim_ctrl_tbl$Pop_scenario.string <- tolower(
    Sim_ctrl_tbl$Pop_scenario.string
  )
  Sim_ctrl_tbl$Scenario_ID.string <- tolower(
    Sim_ctrl_tbl$Scenario_ID.string
  )

  Start_date <- min(Sim_ctrl_tbl$Scenario_start.real)
  End_date <- max(Sim_ctrl_tbl$Scenario_end.real)
  Sim_time_steps <- seq(
    Start_date,
    End_date,
    by = Sim_ctrl_tbl$Step_length.real[1]
  )

  cat("Loading LULC aggregation scheme from:", LULC_agg_path, "\n")
  Aggregation_scheme <- read_excel(LULC_agg_path)
  Aggregation_scheme$Class_abbreviation <- tolower(
    Aggregation_scheme$Class_abbreviation
  )

  list(
    Sim_ctrl_tbl = Sim_ctrl_tbl,
    Aggregation_scheme = Aggregation_scheme,
    Config_IDs = unique(Sim_ctrl_tbl$Simulation_num.),
    Sim_time_steps = Sim_time_steps,
    Start_date = Start_date,
    End_date = End_date
  )
}

#' Generate per-configuration and per-timestep output paths.
#' @param Config_IDs A vector of unique configuration IDs to generate paths for.
#' @param Sim_ctrl_tbl A data frame representing the simulation control table with scenario details.
#' @param Sim_time_steps A vector of simulation time steps to generate paths for.
#' @param lulcc_input_dir The base directory where the input LULC rasters are stored, organized by configuration ID and time step.
#' @param base_dir The base directory for output files.
#' @param image_dir The subdirectory under base_dir where output PNG images will be saved.
#' @param raster_dir The subdirectory under base_dir where output TIFF rasters will be saved.
#' @param perc_area_data_dir The subdirectory under base_dir where percentage area data JSON files will be saved.
#' @param mask_names Vector of mask names to generate masked output paths for.
#' @return Data frame of input/output paths plus scenario metadata and existence flags.
generate_file_paths <- function(
  Config_IDs,
  Sim_ctrl_tbl,
  Sim_time_steps,
  lulcc_input_dir,
  base_dir,
  image_dir,
  raster_dir,
  perc_area_data_dir,
  mask_names
) {
  lulc_paths <- lapply(Config_IDs, function(Config_ID) {
    Config_details <- Sim_ctrl_tbl[Sim_ctrl_tbl$Simulation_num. == Config_ID, ]

    lulc_time_paths <- as.data.frame(rbindlist(lapply(
      Sim_time_steps,
      function(Time_step) {
        input_path <- file.path(
          lulcc_input_dir,
          Config_ID,
          paste0(
            "simulated_LULC_simID_",
            Config_ID,
            "_year_",
            Time_step,
            ".tif"
          )
        )

        base_name <- paste0(
          "lulc-",
          Time_step,
          "-",
          Config_details$Climate_scenario.string,
          "-",
          Config_details$Econ_scenario.string,
          "-",
          Config_details$Pop_scenario.string,
          "-",
          Config_details$Scenario_ID.string,
          "-full"
        )

        tif_path <- file.path(base_dir, raster_dir, paste0(base_name, ".tif"))
        png_path <- file.path(base_dir, image_dir, paste0(base_name, ".png"))
        perc_area_path <- file.path(
          base_dir,
          perc_area_data_dir,
          paste0(base_name, "-perc_area.json")
        )

        mask_paths_tif <- sapply(mask_names, function(mask_name) {
          gsub("full", mask_name, tif_path, fixed = TRUE)
        })
        mask_paths_png <- sapply(mask_names, function(mask_name) {
          gsub("full", mask_name, png_path, fixed = TRUE)
        })
        mask_paths_perc_area <- sapply(mask_names, function(mask_name) {
          gsub("full", mask_name, perc_area_path, fixed = TRUE)
        })

        time_step_paths <- list(
          input_path,
          tif_path,
          png_path,
          perc_area_path,
          mask_paths_tif,
          mask_paths_png,
          mask_paths_perc_area
        )
        names(time_step_paths) <- c(
          "Path",
          "lulc_path_tif",
          "lulc_path_png",
          "lulc_path_perc_area",
          paste0("lulc_path_tif_", mask_names),
          paste0("lulc_path_png_", mask_names),
          paste0("lulc_path_perc_area_", mask_names)
        )

        time_step_paths
      }
    )))

    lulc_time_paths$Config_ID <- Config_ID
    lulc_time_paths$Time_step <- Sim_time_steps
    lulc_time_paths$Scenario_ID <- Config_details$Scenario_ID.string
    lulc_time_paths$Climate_scenario <- Config_details$Climate_scenario.string
    lulc_time_paths$Econ_scenario <- Config_details$Econ_scenario.string
    lulc_time_paths$Pop_scenario <- Config_details$Pop_scenario.string
    lulc_time_paths$Exists <- file.exists(lulc_time_paths$Path)

    lulc_time_paths
  })

  do.call(rbind, lulc_paths)
}

#' Load a LULC palette mapping from JSON.
#' @param palette_path Path to a JSON file with `raster_value`, `lulc_name`,
#'   and `colour` fields.
#' @return Data frame with one row per class.
load_lulc_palette <- function(palette_path) {
  palette_data <- jsonlite::fromJSON(palette_path)
  required_cols <- c("raster_value", "lulc_name", "colour")

  if (!all(required_cols %in% colnames(palette_data))) {
    stop(
      paste(
        "Palette JSON must contain columns:",
        paste(required_cols, collapse = ", ")
      )
    )
  }

  palette_data$raster_value <- as.numeric(palette_data$raster_value)
  palette_data
}

#' Create the LULC raster attribute table and colour lookup.
#' @param lulc_raster_path The file path to the LULC raster for which the RAT is to be created.
#' @param Aggregation_scheme A data frame containing the mapping of raster values to class names, with columns 'Aggregated_ID' and 'Aggregated_class_short'.
#' @param palette_data Data frame with `raster_value`, `lulc_name`, and `colour`.
#' @return List with `LULC_rat`, `col_map`, `mask_values`, `glacier_val`,
#'   and `rock_val`.
create_lulc_rat <- function(
  lulc_raster_path,
  Aggregation_scheme,
  palette_data
) {
  raster_ids <- sort(unique(values(rast(lulc_raster_path))))
  raster_ids <- raster_ids[!is.na(raster_ids)]
  rat_ids <- sort(unique(c(raster_ids, 20, 21)))

  LULC_rat <- data.frame(ID = rat_ids)
  LULC_rat$lulc_name <- sapply(LULC_rat$ID, function(id) {
    if (id == 20) {
      return("Lake")
    }
    if (id == 21) {
      return("River")
    }

    matched_class <- unique(unlist(Aggregation_scheme[
      Aggregation_scheme$Aggregated_ID == id,
      "Aggregated_class_short"
    ]))
    if (length(matched_class) == 0) {
      return(paste0("Unknown_", id))
    }

    matched_class[1]
  })

  LULC_rat$colour <- sapply(seq_len(nrow(LULC_rat)), function(i) {
    matched_colour <- palette_data$colour[
      palette_data$raster_value == LULC_rat$ID[i]
    ]

    if (length(matched_colour) == 0) {
      matched_colour <- palette_data$colour[
        palette_data$lulc_name == as.character(LULC_rat$lulc_name[i])
      ]
    }

    if (length(matched_colour) == 0) {
      stop(
        paste(
          "No palette entry found for raster value",
          LULC_rat$ID[i],
          "and class",
          as.character(LULC_rat$lulc_name[i])
        )
      )
    }

    matched_colour[1]
  })
  LULC_rat <- LULC_rat[order(LULC_rat$ID), ]
  LULC_rat$lulc_name <- factor(
    LULC_rat$lulc_name,
    levels = LULC_rat$lulc_name,
    ordered = TRUE
  )

  mask_values <- unlist(Aggregation_scheme[
    Aggregation_scheme$NOAS04_class_ENG %in%
      c("Lakes", "Rivers", "Glaciers, perpetual snow"),
    "NOAS04_ID"
  ])
  names(mask_values) <- c(20, 21, 19)

  list(
    LULC_rat = LULC_rat,
    col_map = setNames(LULC_rat$colour, as.character(LULC_rat$ID)),
    mask_values = mask_values,
    glacier_val = unlist(Aggregation_scheme[
      Aggregation_scheme$NOAS04_class_ENG == "Glaciers, perpetual snow",
      "NOAS04_ID"
    ]),
    rock_val = unlist(Aggregation_scheme[
      Aggregation_scheme$NOAS04_class_ENG == "Rocks",
      "NOAS04_ID"
    ])
  )
}

#' Load the reference LULC raster and flatten it to a data frame.
#' @param Non_agg_lulc_path The file path to the non-aggregated LULC raster to be loaded.
#' @param ProjCH The target coordinate reference system (CRS) to which the LULC raster should be projected if it is not already in that CRS.
#' @return List with the raster, a coordinate-enriched data frame, and the name
#'   of the value column.
load_reference_lulc <- function(Non_agg_lulc_path, ProjCH) {
  ref_LULC <- rast(Non_agg_lulc_path)
  if (crs(ref_LULC) != ProjCH) {
    ref_LULC <- project(ref_LULC, ProjCH)
  }

  ref_LULC_dat <- terra::as.data.frame(ref_LULC, na.rm = FALSE)
  ref_LULC_dat$ID <- seq.int(nrow(ref_LULC))
  ref_LULC_dat <- cbind(ref_LULC_dat, crds(ref_LULC, na.rm = FALSE))

  list(
    ref_LULC = ref_LULC,
    ref_LULC_dat = ref_LULC_dat,
    value_col = setdiff(colnames(ref_LULC_dat), c("ID", "x", "y"))[1]
  )
}

#' Load glacier scenario IDs for a given scenario and start date.
#' @param glacier_scenario_dir The directory path where the glacier scenario files are stored.
#' @param Scenario A string representing the scenario identifier to match in the glacier scenario files.
#' @param Start_date Start date used to select the glacier index column.
#' @return List with `Glacier_IDs` and `Non_glacier_IDs`.
load_glacier_scenario <- function(glacier_scenario_dir, Scenario, Start_date) {
  glacier_files <- list.files(glacier_scenario_dir, full.names = TRUE)
  scenario_glacier_file <- glacier_files[grepl(
    Scenario,
    glacier_files,
    ignore.case = TRUE
  )]

  if (length(scenario_glacier_file) != 1) {
    stop(
      paste(
        "Expected exactly one glacier file for scenario",
        Scenario,
        "but found",
        length(scenario_glacier_file)
      )
    )
  }

  glacier_index <- readRDS(scenario_glacier_file)[, c(
    "ID_loc",
    paste(Start_date)
  )]

  list(
    Glacier_IDs = glacier_index[
      glacier_index[[paste(Start_date)]] == 1,
      "ID_loc"
    ],
    Non_glacier_IDs = glacier_index[
      glacier_index[[paste(Start_date)]] == 0,
      "ID_loc"
    ]
  )
}

#' Create a scenario-specific LULC raster from glacier status.
#' @param ref_LULC_dat Data frame version of the reference raster.
#' @param value_col Name of the LULC value column in `ref_LULC_dat`.
#' @param Glacier_IDs A vector of IDs corresponding to locations classified as glaciers in the
#' glacier scenario data.
#' @param Non_glacier_IDs A vector of IDs corresponding to locations classified as non glaciers in the glacier scenario data.
#' @param glacier_val The numeric value to assign to locations identified as glaciers in the scenario.
#' @param rock_val The numeric value to assign to locations identified as non-glaciers in the scenario.
#' @param ProjCH Coordinate reference system for the returned raster.
#' @return Scenario-specific `SpatRaster`.
create_scenario_lulc <- function(
  ref_LULC_dat,
  value_col,
  Glacier_IDs,
  Non_glacier_IDs,
  glacier_val,
  rock_val,
  ProjCH
) {
  scenario_LULC_dat <- ref_LULC_dat

  scenario_LULC_dat[
    scenario_LULC_dat$ID %in% Non_glacier_IDs,
    value_col
  ] <- rock_val
  scenario_LULC_dat[
    scenario_LULC_dat$ID %in% Glacier_IDs,
    value_col
  ] <- glacier_val
  scenario_LULC_dat[
    which(
      scenario_LULC_dat[[value_col]] == glacier_val &
        !(scenario_LULC_dat$ID %in% Glacier_IDs)
    ),
    value_col
  ] <- rock_val

  rast(scenario_LULC_dat[, c("x", "y", value_col)], crs = ProjCH)
}

#' Apply scenario masking to a single timestep raster.
#' @param lulc_path The file path to the LULC raster for the current time step to be processed.
#' @param scenario_LULC A SpatRaster object representing the scenario LULC raster that will be used for masking the LULC layer.
#' @param time_step The current time step being processed.
#' @param first_time_step The first time step in the simulation, used to determine which mask values to apply.
#' @param mask_values Named vector of mask values to update from `scenario_LULC`.
#' @param ProjCH Coordinate reference system for the timestep raster.
#' @return Processed `SpatRaster`.
process_timestep <- function(
  lulc_path,
  scenario_LULC,
  time_step,
  first_time_step,
  mask_values,
  ProjCH
) {
  lulc_layer <- rast(lulc_path)
  crs(lulc_layer) <- ProjCH

  if (time_step == first_time_step) {
    time_step_mask_values <- mask_values
  } else {
    time_step_mask_values <- mask_values[names(mask_values) %in% c(20, 21)]
  }

  for (j in seq_along(time_step_mask_values)) {
    lulc_layer <- terra::mask(
      x = lulc_layer,
      mask = scenario_LULC,
      maskvalues = time_step_mask_values[j],
      updatevalue = as.numeric(names(time_step_mask_values)[j])
    )
  }

  lulc_layer
}

#' Load a mask layer from a shapefile or raster.
#' @param mask_path The file path to the mask layer to be loaded, which can be either a shapefile (with a .shp extension) or a GeoTIFF raster (with a .tif extension).
#' @return A `SpatVector` or `SpatRaster`.
load_mask_layer <- function(mask_path) {
  if (!file.exists(mask_path)) {
    stop(paste("Mask file does not exist:", mask_path))
  }
  if (grepl("\\.shp$", mask_path, ignore.case = TRUE)) {
    return(vect(mask_path))
  }
  if (grepl("\\.tif$", mask_path, ignore.case = TRUE)) {
    return(rast(mask_path))
  }
  stop(paste("Unsupported mask file type for:", mask_path))
}

#' Build area summaries for a main raster and optional focus regions.
#' @param main_raster A SpatRaster object representing the main LULC raster for which area results are to be computed.
#' @param main_region_name A string representing the name of the main region (e.g., "full extent") to be used as a label in the area results.
#' @param Aggregation_scheme Data frame used to translate raster IDs to class names.
#' @param focus_regions Optional named list of additional raster masks.
#' @return Named list of area summaries.
build_area_results <- function(
  main_raster,
  main_region_name,
  Aggregation_scheme,
  focus_regions = NULL
) {
  area_results <- list()
  area_results[[main_region_name]] <- process_lulc_area(
    raster = main_raster,
    Aggregation_scheme = Aggregation_scheme,
    context_label = main_region_name
  )

  if (!is.null(focus_regions) && length(focus_regions) > 0) {
    for (area_name in names(focus_regions)) {
      focus_mask_path <- focus_regions[[area_name]]
      if (!file.exists(focus_mask_path)) {
        cat(
          "Warning: focus mask file not found for area '",
          area_name,
          "': ",
          focus_mask_path,
          " - skipping.\n",
          sep = ""
        )
        next
      }

      focus_raster <- terra::mask(main_raster, rast(focus_mask_path))
      area_results[[area_name]] <- process_lulc_area(
        raster = focus_raster,
        Aggregation_scheme = Aggregation_scheme,
        context_label = area_name
      )
    }
  }

  Filter(Negate(is.null), area_results)
}

#' Write a list to a pretty-printed JSON file.
#' @param data_list Named list to serialize.
#' @param output_path Output JSON path.
#' @return No return value. Writes a file to disk.
save_area_json <- function(data_list, output_path) {
  write(toJSON(data_list, pretty = TRUE), file = output_path)
}

#' Save raster, area summary JSON, and PNG outputs for one map.
#' @param raster_layer A SpatRaster object representing the LULC raster layer to be saved as a TIFF file and used for generating the indexed PNG image.
#' @param area_results A list containing the area results to be saved as a JSON file at the specified area path.
#' @param tif_path The file path where the LULC raster layer should be saved as a TIFF file.
#' @param png_path The file path where the indexed PNG image should be saved.
#' @param area_path The file path where the area results should be saved as a JSON file.
#' @param col_map A named vector mapping raster values (as character) to colour codes, used for generating the indexed PNG image.
#' @param overwrite A boolean indicating whether to overwrite existing files.
#' @return No return value.
save_lulc_outputs <- function(
  raster_layer,
  area_results,
  tif_path,
  png_path,
  area_path,
  col_map,
  overwrite = TRUE
) {
  writeRaster(raster_layer, filename = tif_path, overwrite = overwrite)
  save_area_json(area_results, area_path)
  save_indexed_png(
    raster_obj = raster_layer,
    output_path = png_path,
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
    verbose = FALSE
  )
}

#' Calculate percentage area changes between two saved area summaries.
#' @param current_perc_area_path Path to the current percentage-area JSON.
#' @param first_perc_area_path Path to the baseline percentage-area JSON.
#' @return Named list of regional class-wise percentage changes.
calculate_area_changes <- function(
  current_perc_area_path,
  first_perc_area_path
) {
  first_all_regions <- fromJSON(first_perc_area_path, simplifyDataFrame = FALSE)
  current_all_regions <- fromJSON(
    current_perc_area_path,
    simplifyDataFrame = FALSE
  )

  all_region_changes <- list()

  for (region_name in names(first_all_regions)) {
    if (!region_name %in% names(current_all_regions)) {
      cat(
        "Warning: Region",
        region_name,
        "not found in current data, skipping\n"
      )
      next
    }

    first_region_data <- first_all_regions[[region_name]]
    current_region_data <- current_all_regions[[region_name]]

    first_data <- data.frame(
      class_name = names(first_region_data),
      perc_area = as.numeric(unlist(first_region_data)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    current_data <- data.frame(
      class_name = names(current_region_data),
      perc_area = as.numeric(unlist(current_region_data)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    area_change <- merge(
      first_data,
      current_data,
      by = "class_name",
      suffixes = c("_start", "_end")
    )
    area_change$perc_area_change <- ((area_change$perc_area_end -
      area_change$perc_area_start) /
      area_change$perc_area_start) *
      100
    area_change <- area_change[
      order(match(area_change$class_name, first_data$class_name)),
    ]

    all_region_changes[[region_name]] <- setNames(
      as.list(area_change$perc_area_change),
      area_change$class_name
    )
  }

  all_region_changes
}

#' Process and save outputs for one masked extent.
#' @param lulc_layer Processed timestep raster.
#' @param scenario_row One-row data frame of output paths and metadata.
#' @param mask_name Name of the current mask.
#' @param mask_path Path to the mask layer.
#' @param Aggregation_scheme Aggregation table used for area summaries.
#' @param focus_regions Optional named list of extra focus-region masks.
#' @param col_map Named raster-value to colour lookup.
#' @param overwrite Whether to overwrite existing outputs.
#' @return No return value.
process_masked_output <- function(
  lulc_layer,
  scenario_row,
  mask_name,
  mask_path,
  Aggregation_scheme,
  focus_regions,
  col_map,
  overwrite = TRUE
) {
  mask_path_tif <- scenario_row[[paste0("lulc_path_tif_", mask_name)]]
  map_path_png <- scenario_row[[paste0("lulc_path_png_", mask_name)]]
  mask_path_data <- scenario_row[[paste0("lulc_path_perc_area_", mask_name)]]

  if (
    !overwrite &&
      file.exists(mask_path_tif) &&
      file.exists(map_path_png) &&
      file.exists(mask_path_data)
  ) {
    message(paste(
      "Files already exist, skipping:",
      mask_path_tif,
      map_path_png,
      mask_path_data
    ))
    return(invisible(NULL))
  }

  mask_layer <- load_mask_layer(mask_path)
  cropped_lulc <- terra::crop(x = lulc_layer, y = mask_layer)
  masked_lulc <- terra::mask(
    x = cropped_lulc,
    mask = mask_layer,
    updatevalue = NA
  )

  area_results <- build_area_results(
    main_raster = masked_lulc,
    main_region_name = mask_name,
    Aggregation_scheme = Aggregation_scheme,
    focus_regions = focus_regions
  )

  save_lulc_outputs(
    raster_layer = masked_lulc,
    area_results = area_results,
    tif_path = mask_path_tif,
    png_path = map_path_png,
    area_path = mask_path_data,
    col_map = col_map,
    overwrite = TRUE
  )
}

#' Process and save outputs for the full unmasked extent.
#' @param lulc_layer Processed timestep raster.
#' @param scenario_row One-row data frame of output paths and metadata.
#' @param Aggregation_scheme Aggregation table used for area summaries.
#' @param focus_regions Optional named list of extra focus-region masks.
#' @param col_map Named raster-value to colour lookup.
#' @param overwrite Whether to overwrite existing outputs.
#' @return No return value.
process_full_extent_output <- function(
  lulc_layer,
  scenario_row,
  Aggregation_scheme,
  focus_regions,
  col_map,
  overwrite = TRUE
) {
  if (
    !overwrite &&
      file.exists(scenario_row$lulc_path_tif) &&
      file.exists(scenario_row$lulc_path_png) &&
      file.exists(scenario_row$lulc_path_perc_area)
  ) {
    message(paste(
      "Files already exist, skipping:",
      scenario_row$lulc_path_tif,
      scenario_row$lulc_path_png,
      scenario_row$lulc_path_perc_area
    ))
    return(invisible(NULL))
  }

  area_results <- build_area_results(
    main_raster = lulc_layer,
    main_region_name = "full",
    Aggregation_scheme = Aggregation_scheme,
    focus_regions = focus_regions
  )

  save_lulc_outputs(
    raster_layer = lulc_layer,
    area_results = area_results,
    tif_path = scenario_row$lulc_path_tif,
    png_path = scenario_row$lulc_path_png,
    area_path = scenario_row$lulc_path_perc_area,
    col_map = col_map,
    overwrite = TRUE
  )
}

#' Calculate and save area-change JSON files for one configuration.
#' @param config_lulc_df Data frame of timestep paths for one configuration.
#' @param mask_names Character vector of mask names to process.
#' @param base_dir Base output directory.
#' @param area_chg_data_dir Relative output directory for area-change JSON.
#' @param overwrite Whether to overwrite existing outputs.
#' @return No return value.
calculate_and_save_area_changes <- function(
  config_lulc_df,
  mask_names,
  base_dir,
  area_chg_data_dir,
  overwrite = TRUE
) {
  first_time_step <- min(config_lulc_df$Time_step)
  config_time_steps <- unique(config_lulc_df$Time_step)
  config_time_steps <- config_time_steps[config_time_steps != first_time_step]

  for (time_step in config_time_steps) {
    for (mask_name in mask_names) {
      area_change_path <- file.path(
        base_dir,
        area_chg_data_dir,
        paste0(
          "lulc-",
          time_step,
          "-",
          config_lulc_df$Climate_scenario[
            config_lulc_df$Time_step == time_step
          ],
          "-",
          config_lulc_df$Econ_scenario[config_lulc_df$Time_step == time_step],
          "-",
          config_lulc_df$Pop_scenario[config_lulc_df$Time_step == time_step],
          "-",
          config_lulc_df$Scenario_ID[config_lulc_df$Time_step == time_step],
          "-",
          mask_name,
          "-perc_area_chg.json"
        )
      )

      if (!overwrite && file.exists(area_change_path)) {
        message(paste(
          "Area change file already exists, skipping:",
          area_change_path
        ))
        next
      }

      perc_area_col <- paste0("lulc_path_perc_area_", mask_name)
      if (!perc_area_col %in% colnames(config_lulc_df)) {
        stop(paste("No column for", perc_area_col, "in config_lulc_df"))
      }

      first_mask_path <- config_lulc_df[
        config_lulc_df$Time_step == first_time_step,
        perc_area_col
      ]
      current_mask_path <- config_lulc_df[
        config_lulc_df$Time_step == time_step,
        perc_area_col
      ]

      area_changes <- calculate_area_changes(
        current_perc_area_path = current_mask_path,
        first_perc_area_path = first_mask_path
      )
      save_area_json(area_changes, area_change_path)
    }
  }
}

#' Prepare LULC rasters, PNGs, and chart JSON for all scenarios.
#' @param lulcc_input_dir The directory path where the input LULC files are stored.
#' @param image_dir Output subdirectory for PNGs.
#' @param raster_dir Output subdirectory for rasters.
#' @param perc_area_data_dir Output subdirectory for area-summary JSON.
#' @param area_chg_data_dir Output subdirectory for area-change JSON.
#' @param base_dir Base output directory.
#' @param Sim_ctrl_tbl_path Path to the simulation control table.
#' @param ProjCH Coordinate reference system for raster processing.
#' @param LULC_agg_path Path to the LULC aggregation scheme.
#' @param palette_data Data frame with palette metadata loaded from JSON.
#' @param Non_agg_lulc_path Path to the non-aggregated reference LULC raster.
#' @param glacier_scenario_dir Directory containing glacier scenario indices.
#' @param Use_parallel Whether to process scenarios in parallel.
#' @param num_workers Number of workers for parallel processing.
#' @param map_masks Optional named list of masks to apply for output maps.
#' @param focus_regions Optional named list of masks used only for area summaries.
#' @param overwrite Whether to overwrite existing outputs.
#' @return NULL
prepare_lulc_files <- function(
  lulcc_input_dir,
  image_dir,
  raster_dir,
  perc_area_data_dir,
  area_chg_data_dir,
  base_dir,
  Sim_ctrl_tbl_path,
  ProjCH,
  LULC_agg_path,
  palette_data,
  Non_agg_lulc_path,
  glacier_scenario_dir,
  Use_parallel,
  num_workers,
  map_masks,
  focus_regions,
  overwrite = TRUE
) {
  ensure_dir(file.path(base_dir, image_dir))
  ensure_dir(file.path(base_dir, raster_dir))
  ensure_dir(file.path(base_dir, perc_area_data_dir))
  ensure_dir(file.path(base_dir, area_chg_data_dir))

  config_data <- load_and_prepare_config(
    Sim_ctrl_tbl_path = Sim_ctrl_tbl_path,
    LULC_agg_path = LULC_agg_path
  )

  mask_names <- if (is.null(map_masks)) "full" else names(map_masks)

  cat(
    "Building file paths for",
    length(config_data$Config_IDs),
    "configurations\n"
  )
  lulc_df <- generate_file_paths(
    Config_IDs = config_data$Config_IDs,
    Sim_ctrl_tbl = config_data$Sim_ctrl_tbl,
    Sim_time_steps = config_data$Sim_time_steps,
    lulcc_input_dir = lulcc_input_dir,
    base_dir = base_dir,
    image_dir = image_dir,
    raster_dir = raster_dir,
    perc_area_data_dir = perc_area_data_dir,
    mask_names = mask_names
  )

  if (!all(lulc_df$Exists)) {
    cat("Missing files:\n")
    print(lulc_df$Path[!lulc_df$Exists])
    stop("Some LULC files do not exist. Please check the paths.")
  }

  lulc_rat_data <- create_lulc_rat(
    lulc_raster_path = lulc_df$Path[1],
    Aggregation_scheme = config_data$Aggregation_scheme,
    palette_data = palette_data
  )
  reference_lulc_data <- load_reference_lulc(
    Non_agg_lulc_path = Non_agg_lulc_path,
    ProjCH = ProjCH
  )

  if (Use_parallel) {
    plan(multisession, workers = num_workers)
  } else {
    plan(sequential)
  }

  future_sapply(unique(lulc_df$Scenario_ID), function(Scenario) {
    cat("\n### Processing scenario:", Scenario, "###\n")

    glacier_data <- load_glacier_scenario(
      glacier_scenario_dir = glacier_scenario_dir,
      Scenario = Scenario,
      Start_date = config_data$Start_date
    )

    scenario_LULC <- create_scenario_lulc(
      ref_LULC_dat = reference_lulc_data$ref_LULC_dat,
      value_col = reference_lulc_data$value_col,
      Glacier_IDs = glacier_data$Glacier_IDs,
      Non_glacier_IDs = glacier_data$Non_glacier_IDs,
      glacier_val = lulc_rat_data$glacier_val,
      rock_val = lulc_rat_data$rock_val,
      ProjCH = ProjCH
    )

    scenario_lulc_df <- lulc_df[lulc_df$Scenario_ID == Scenario, ]
    first_time_step <- min(scenario_lulc_df$Time_step)

    for (i in seq_len(nrow(scenario_lulc_df))) {
      scenario_row <- scenario_lulc_df[i, ]

      if (!scenario_row$Exists) {
        message(paste("File does not exist:", scenario_row$Path))
        next
      }

      lulc_layer <- process_timestep(
        lulc_path = scenario_row$Path,
        scenario_LULC = scenario_LULC,
        time_step = scenario_row$Time_step,
        first_time_step = first_time_step,
        mask_values = lulc_rat_data$mask_values,
        ProjCH = ProjCH
      )

      if (!is.null(map_masks) && length(map_masks) > 0) {
        for (mask_name in mask_names) {
          process_masked_output(
            lulc_layer = lulc_layer,
            scenario_row = scenario_row,
            mask_name = mask_name,
            mask_path = map_masks[[mask_name]],
            Aggregation_scheme = config_data$Aggregation_scheme,
            focus_regions = focus_regions,
            col_map = lulc_rat_data$col_map,
            overwrite = overwrite
          )
        }
      } else {
        process_full_extent_output(
          lulc_layer = lulc_layer,
          scenario_row = scenario_row,
          Aggregation_scheme = config_data$Aggregation_scheme,
          focus_regions = focus_regions,
          col_map = lulc_rat_data$col_map,
          overwrite = overwrite
        )
      }
    }

    for (Config_ID in unique(scenario_lulc_df$Config_ID)) {
      config_lulc_df <- scenario_lulc_df[
        scenario_lulc_df$Config_ID == Config_ID,
      ]

      calculate_and_save_area_changes(
        config_lulc_df = config_lulc_df,
        mask_names = mask_names,
        base_dir = base_dir,
        area_chg_data_dir = area_chg_data_dir,
        overwrite = overwrite
      )
    }
  })

  cat("\n=== All scenarios processed successfully ===\n")
}

### =========================================================================
### Main
### =========================================================================

cat("\n=== Starting LULC Data Preparation Script ===\n")
cat("Loading configuration file\n")
config_file <- Sys.getenv("FUTURE_EI_CONFIG_FILE")
cat("Config file path:", config_file, "\n")
if (!file.exists(config_file)) {
  stop(paste0(
    "Config file FUTURE_EI_CONFIG_FILE (",
    config_file,
    ") does not exist."
  ))
}
config <- yaml.load_file(config_file)
cat("Configuration loaded successfully\n")
bash_vars <- config$bash_variables
config <- config$Summarisation
cat("Using summarisation configuration\n")

cat("\n=== Setting up directories ===\n")
if (!is.null(config$OutputDir) && nchar(config$OutputDir) > 0) {
  web_platform_dir <- config$OutputDir
  cat("Using OutputDir from config:", web_platform_dir, "\n")
} else {
  web_platform_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform"
  cat("Using default web platform directory:", web_platform_dir, "\n")
}

ensure_dir(web_platform_dir)

Mask_dir <- file.path(web_platform_dir, "Masks")
ensure_dir(Mask_dir)

palette_path <- config$LULCPalettePath
if (!grepl("^(/|[A-Za-z]:)", palette_path)) {
  palette_path <- file.path(dirname(config_file), palette_path)
}
LULC_palette <- load_lulc_palette(palette_path)

cat("\n=== Starting LULC file preparation ===\n")
options(future.globals.maxSize = 8 * 1024^3)

focus_regions <- NULL
if (
  !is.null(config$PercAreaFocusRegions) &&
    length(config$PercAreaFocusRegions) > 0
) {
  cat("\n=== Loading focus regions from config ===\n")
  focus_regions <- config$PercAreaFocusRegions

  for (region_name in names(focus_regions)) {
    region_path <- focus_regions[[region_name]]
    if (!grepl("^(/|[A-Za-z]:)", region_path)) {
      focus_regions[[region_name]] <- file.path(web_platform_dir, region_path)
    }
    cat("  -", region_name, ":", focus_regions[[region_name]], "\n")
  }
} else {
  cat("\n=== No focus regions specified in config ===\n")
}

prepare_lulc_files(
  lulcc_input_dir = config$LULCCInputDir,
  image_dir = config$ImageDir,
  raster_dir = config$RasterDir,
  perc_area_data_dir = config$PercAreaDir,
  area_chg_data_dir = config$AreaChgDir,
  base_dir = web_platform_dir,
  Sim_ctrl_tbl_path = Sys.getenv("LULCC_M_SIM_CONTROL_TABLE"),
  ProjCH = config$ProjCH,
  LULC_agg_path = config$LULCAggPath,
  palette_data = LULC_palette,
  Non_agg_lulc_path = config$LULCRefRaster,
  glacier_scenario_dir = config$GlacierScenarioDir,
  Use_parallel = config$Parallel,
  num_workers = config$NWorkers,
  map_masks = config$Masks,
  focus_regions = focus_regions,
  overwrite = TRUE
)
