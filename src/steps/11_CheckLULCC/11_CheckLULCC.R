#' Perform intensity analysis on simulated LULC layers
#'
#' Input: Folder name of simulation maps
#' for the set of maps for each simulation:
#'     Perform Interval level intensity anslysis
#'     save results as rds
#' Output: Intensity analysis results 
#'
#' @environment focal_lulc
#' @config $FUTURE_EI_CONFIG_FILE (yaml file)
#' @date 2024-07-16,
#' @author Benjamin Black, Carlson Büth
#' @note This script uses code by OpenLand
#' @references Reginal Exavier and Peter Zeilhofer. OpenLand: Software for
#' Quantitative Analysis and Visualization of Land Use and Cover Change.
#' The R Journal, v. 12, n. 2, p. 359–371, 2021.
#' https://doi.org/10.32614/RJ-2021-021.
#'
#'
#' @docType script
#'


# Load libraries
packs <- c("raster", "yaml", "future", "future.apply", "stringr", "readxl")
invisible(lapply(packs, require, character.only = TRUE))

options(future.rng.onMisuse = "ignore")

## Functions

#' acc_changes: 
#' 
#' Count cumulative pixel wise changes in an LULC raster time series
#'
#' This function calculates the number of times a pixel undergoes changes during
#' a time series. It returns a raster with the number of changes as
#' pixel value and a table containing the areal percentage of every pixel value
#' (number of changes).
#'
#'
#' @param input_rasters The path for the Raster* directory or list of Raster* to
#' be analysed.
#'
#'
#' @return DF: frequency and % of cumulative pixel-wise change values
#' @export
#'
#'
#'
#' @examples
#' \donttest{
#' url <- "https://zenodo.org/record/3685230/files/SaoLourencoBasin.rda?download=1"
#' temp <- tempfile()
#' download.file(url, temp, mode = "wb")
#' # downloading the SaoLourencoBasin dataset
#' load(temp)
#' # the acc_changes() function, with the SaoLourencoBasin dataset
#' acc_changes(SaoLourencoBasin)
#' }
#'

acc_changes <- function(input_rasters) {

  # Convert rasterstack to list
  rList <- raster::unstack(input_rasters)

  # get number of rasters in list
  n_raster <- length(rList)

  # Stop if insufficient number of rasters
  if (n_raster < 2) {
    stop('acc_changes needs at least 2 rasters')
  }

  #calculate pairwise differences between rasters
  difflist <- mapply(
    function(x, y)
      raster::overlay(
        x,
        y,
        fun = function(x1, x2)
          ifelse((x1 != x2), 1, 0)
      ),
    x = rList[1:(length(rList) - 1)],
    y = rList[2:length(rList)],
    SIMPLIFY = FALSE
  )

  #sum rasters in stack
  sumraster <- sum(raster::stack(difflist))

  # empty values for columns
  Freq <- Var1 <- NULL

  # summarize frequencies of changes
  df01_values <- table(matrix(sumraster))

  # change column types and names
  df_values <- dplyr::mutate(data.frame(df01_values),
                             Var1 = as.character(Var1),
                             Var1 = as.integer(Var1),
                             Percent = Freq / sum(Freq) * 100)
  names(df_values) <- c("Cumulative_changes", "Npixels", "Percent")
  return(df_values)

}

#' intensityAnalysis:
#' 
#' Calculate contingency tables over a time series of rasters and then
#' calculate measures of intensity analysis
#'
#' @param input_rasters Rasterstack* object. See \cr
#' \code{\link[raster]{raster}} for more information about supported file types.
#' @param pixelresolution numeric. The pixel spatial resolution in meter.
#'
#' @docType methods
#' @import dplyr
#'
#' @return A list that contains 5 objects.
#' \itemize{
#'   \item \code{lulc_Mulstistep}: \code{<tibble>} Contingency table for all
#'   analysed time steps, containing 8 columns:
#'   \enumerate{
#'   \item Period: \code{<chr>} The period \emph{[Yt, Yt+1]}.
#'   \item From: \code{<dbl>} numerical code of a LUC category \emph{i}.
#'   \item To: \code{<dbl>} numerical code of a LUC category \emph{j}.
#'   \item km2: \code{<dbl>} Area in square kilometers that transited from the
#'   category \emph{i}
#'    to category \emph{j} in the period from \emph{Yt} to \emph{Yt+1}.
#'   \item Interval: \code{<dbl>} Interval of years between the first and
#'    the last year of the period \emph{[Yt, Yt+1]}.
#'   \item QtPixel: \code{<int>} Pixel count that transited from the categories
#'    \emph{i}
#'    to category \emph{j} in the period from \emph{Yt} to \emph{Yt+1}.
#'   \item yearFrom: \code{<chr>} The year that the change comes from \emph{[Yt]}.
#'   \item yearTo: \code{<chr>} The year that the change goes for \emph{[Yt+1]}.
#'   }
#'   \item \code{lulc_Onestep}:\code{<tibble>} Contingency table for the entire
#'   analysed period \emph{[Y1, YT]}, containing
#'   8 columns identical with \code{lulc_Mulstistep}.
#'   \item \code{tb_legend}: \code{<tibble>} A table of the pixel value, his
#'   name and color containing 3 columns:
#'   \enumerate{
#'   \item categoryValue: \code{<dbl>} the pixel value of the LUC category.
#'   \item categoryName: \code{<factor>} randomly created string associated with
#'    a given pixel value of a LUC category.
#'   \item color: \code{<chr>} random color associated with the given pixel value
#'    of a LUC category.
#'   Before further analysis, one would like to change the \code{categoryName}
#'   and \code{color} values.
#'     \itemize{
#'         \item Therefore the category names have to be in the same order as the
#'          \code{categoryValue}
#'         and the \code{levels} should be put in the right order for legend
#'         plotting. Like:
#'         \preformatted{
#'
#'         myobject$tb_legend$categoryName <- factor(c("name1", "name2", "name3", "name4"),
#'                                                levels = c("name3", "name2", "name1", "name4"))}
#'         \item The colors have to in the same order as the values in the \code{categoryValue} column. Colors can be given by the
#'        color name (eg. "black") or an HEX value (eg. #FFFFFF). Like:
#'        \preformatted{
#'
#'        myobject$tb_legend$color <- c("#CDB79E", "red", "#66CD00", "yellow")}}}
#'   \item \code{totalArea}: \code{<tibble>}  A table with the total area of the
#'    study area containing 2 columns:
#'   \enumerate{
#'   \item area_km2: \code{<numeric>} The total area in square kilometers.
#'   \item QtPixel: \code{<numeric>} The total area in pixel counts.
#'   }
#'   \item \code{totalInterval}: \code{<numeric>} Total interval of the analysed
#'    time series in years.
#'   }
#'   
#' Performs the intensity analysis based on cross-tabulation matrices of each
#' time step
#'
#' This function implements an Intensity Analysis (IA) according to Aldwaik &
#' Pontius (2012), a quantitative method to analyze time series of land use and
#' cover (LUC) maps.  For IA, a cross-tabulation matrix is composed for each LUC
#' transition step in time.
#'
#' IA includes three levels of analysis of LUC changes. Consecutive analysis
#' levels detail hereby information given by the previous analysis level
#' \cite{(Aldwaik and Pontius, 2012, 2013)}.
#'
#'
#' \enumerate{
#'  \item The \emph{interval level} examines how the size and speed of change
#'  vary across time intervals.
#'  \item The \emph{category level} examines how the size and intensity of gross
#'  losses and gross gains in each category vary across categories for each time
#'  interval.
#'  \item The \emph{transition level} examines how the size and intensity of a
#'  category’s transitions vary across the other categories that are available
#'  for that transition.
#'   }
#'
#'
#' At each analysis level, the method tests for stationarity of patterns across
#' time intervals.
#'
#' \bold{The function returns a list with 6 objects:}
#' \enumerate{
#'  \item lulc_table: \code{tibble}. Contingency table of LUC transitions at all
#'   analysed time steps, containing 6 columns:
#'    \enumerate{
#'    \item Period:  \code{<fct>}. Evaluated period of transition in the format
#'    \code{year t - year t+1}.
#'    \item From: \code{<fct>}. The category in year t.
#'    \item To: \code{<fct>}. The category in year t+1.
#'    \item km2: \code{<dbl>}. Area in square kilometers that transited from the
#'     category \code{From}.
#'    to the category \code{To} in the period.
#'    \item QtPixel: \code{<int>}. Number of pixels that transited from.
#'    the category \code{From} to the category \code{To} in the period.
#'    \item Interval: \code{<int>}. Interval in years of the evaluated period.
#'
#'    }
#'
#'  \item \emph{lv1_tbl}: An \code{\linkS4class{Interval}} object containing the
#'   \emph{St} and \emph{U} values.
#'  \item \emph{category_lvlGain}: A \code{\linkS4class{Category}} object
#'  containing the gain of the LUC category in a period (\emph{Gtj}).
#'  \item \emph{category_lvlLoss}: A \code{\linkS4class{Category}} object
#'  containing the loss of the LUC category in a period (\emph{Lti}).
#'
#'   }
#'
#' @references Aldwaik, S. Z. and Pontius, R. G. (2012) ‘Intensity analysis to unify
#' measurements of size and stationarity of land changes by interval, category, and
#' transition, Landscape and Urban Planning. Elsevier B.V., 106(1), pp. 103–114.
#' \doi{10.1016/j.landurbplan.2012.02.010}.
#'
#' Aldwaik, S. Z. and Pontius, R. G. (2013) ‘Map errors that could account for deviations
#' from a uniform intensity of land change, International Journal of Geographical
#' Information Science. Taylor & Francis, 27(9), pp. 1717–1739. \doi{10.1080/13658816.2013.787618}.
#'
#'
#' @export
#'
#' @importFrom raster unstack crosstab compareRaster raster values stack overlay brick
#'
#' @examples
#' \donttest{
#' url <- "https://zenodo.org/record/3685230/files/SaoLourencoBasin.rda?download=1"
#' temp <- tempfile()
#' download.file(url, temp, mode = "wb") #downloading the online dataset
#' load(temp)
#' # the contingencyTable() function, with the SaoLourencoBasin dataset
#' contingencyTable(input_rasters = SaoLourencoBasin, pixelresolution = 30)
#' }
#'
#'

intensityAnalysis <- function(
  input_rasters,
  pixelresolution = 100
) {

  # helper function to compute the cross table of two rasters and set the column names
  table_cross <- function(x, y) {
    contengency <- raster::crosstab(x, y, long = TRUE, progress = "text")
    contengency %>%
      dplyr::mutate(Year_from = colnames(contengency)[1],
                    Year_to = colnames(contengency)[2]) %>%
      dplyr::rename(
        From = colnames(contengency)[1],
        To = colnames(contengency)[2],
        QtPixel = colnames(contengency)[3]
      ) %>%
      dplyr::mutate(From = as.integer(From), To = as.integer(To))
  }

  # calculate contingency table between initial and final year rasters
  table_one <- table_cross(input_rasters[[1]], input_rasters[[raster::nlayers(input_rasters)]])
  # if there are only two layers then the multi-step table is the same as single step
  if (raster::nlayers(input_rasters) == 2) {
    table_multi <- table_one
  }else { # else if there are > 2 rasters then calculate contingency tables between all pairs of rasters
    input_rasters_multi <- raster::unstack(input_rasters)
    table_multi <- Reduce(rbind,
                          mapply(function(x, y)
                                   table_cross(x, y), input_rasters_multi[1:(length(input_rasters_multi) - 1)],
                                 input_rasters_multi[2:length(input_rasters_multi)], SIMPLIFY = FALSE))
  }

  oneStep <- table_one %>%
    dplyr::arrange(Year_from) %>%
    tidyr::separate(Year_from, c("strings01", "yearFrom"), sep = "_") %>%
    tidyr::separate(Year_to, c("strings02", "yearTo"), sep = "_") %>%
    dplyr::select(-strings01, -strings02) %>%
    dplyr::mutate(yearFrom = as.integer(yearFrom), yearTo = as.integer(yearTo),
                  Interval = yearTo - yearFrom) %>%
    tidyr::unite("Period", c("yearFrom", "yearTo"), sep = "-", remove = FALSE) %>%
    dplyr::select(Period, From, To, QtPixel, Interval, yearFrom, yearTo)

  multiStep <- table_multi %>%
    dplyr::arrange(Year_from) %>%
    tidyr::separate(Year_from, c("strings01", "yearFrom"), sep = "_") %>%
    tidyr::separate(Year_to, c("strings02", "yearTo"), sep = "_") %>%
    dplyr::select(-strings01, -strings02) %>%
    dplyr::mutate(yearFrom = as.integer(yearFrom), yearTo = as.integer(yearTo),
                  Interval = yearTo - yearFrom) %>%
    tidyr::unite("Period", c("yearFrom", "yearTo"), sep = "-", remove = FALSE) %>%
    dplyr::select(Period, From, To, QtPixel, Interval, yearFrom, yearTo)

  # calculate the total time interval and the total pixelValue in the rasters
  allinterval <- as.numeric(dplyr::last(multiStep$yearTo)) - as.numeric(dplyr::first(multiStep$yearFrom))
  areaTotal <- multiStep %>%
    dplyr::group_by(Period) %>%
    dplyr::summarise(QtPixel = sum(QtPixel))

  #convert areaTotal to vector
  totalArea <- unlist(areaTotal[1, "QtPixel"])

  # Select columns
  lulc <- multiStep %>%
    dplyr::select(Period, From, To, QtPixel, Interval)

  # convert period to factor
  lulc$Period <- factor(as.factor(lulc$Period), levels = rev(levels(as.factor(lulc$Period))))

  # ---- Interval analysis (time points) ----
  # EQ1 - St ----

  eq1 <- lulc %>%
    dplyr::filter(From != To) %>%
    dplyr::group_by(Period, Interval) %>%
    dplyr::summarise(intch_QtPixel = sum(QtPixel)) %>% # interval change:intch_QtPixel
    dplyr::mutate(
      PercentChange = (intch_QtPixel / totalArea) * 100,
      St = (intch_QtPixel / (Interval * totalArea)) * 100
    ) %>%
    dplyr::select(1, 4, 5)

  # EQ2 - U ----
  eq2 <- lulc %>%
    dplyr::filter(From != To) %>%
    dplyr::summarise(num02 = sum(QtPixel)) %>% # all area change for the whole period Y1 to YT
    dplyr::mutate(U = (num02 / (allinterval * totalArea)) * 100)

  level01 <- eq1 %>% dplyr::mutate(U = eq2[[2]]) # Type = ifelse(St > U, "Fast", "Slow"))

  # ---- Categorical analysis (LULC class) ----
  # EQ3 - Gtj ----
  num03 <- lulc %>%
    dplyr::filter(From != To) %>%
    dplyr::group_by(Period, To, Interval) %>%
    dplyr::summarise(num03 = sum(QtPixel)) # gross gain category in time point Yt+1

  denom03 <- lulc %>%
    dplyr::group_by(Period, To) %>%
    dplyr::summarise(denom03 = sum(QtPixel)) # total area category in time point Yt+1

  eq3 <- num03 %>%
    dplyr::left_join(denom03, by = c("Period", "To")) %>%
    dplyr::mutate(Gtj = (num03 / (denom03 * Interval)) * 100) %>%
    dplyr::left_join(eq1[c(1, 3)], by = "Period") %>%
    dplyr::select(1, 2, 3, 4, 6, 7) %>%
    dplyr::rename("GG_pixel" = "num03")

  # EQ4 -   Lti ---------
  num04 <- lulc %>%
    dplyr::filter(From != To) %>%
    dplyr::group_by(Period, From, Interval) %>%
    dplyr::summarise(num04 = sum(QtPixel)) # gross loss of category i in time point Yt

  denom04 <- lulc %>%
    dplyr::group_by(Period, From) %>%
    dplyr::summarise(denom04 = sum(QtPixel)) # total area of the category in time point Yt

  eq4 <- num04 %>%
    dplyr::left_join(denom04, by = c("Period", "From")) %>%
    dplyr::mutate(Lti = (num04 / (denom04 * Interval)) * 100) %>%
    dplyr::left_join(eq1[c(1, 3)], by = "Period") %>%
    dplyr::select(1, 2, 3, 4, 6, 7) %>%
    dplyr::rename("GL_pixel" = "num04")

  # ---- Stationarity testing ----

  # ---- Categorical gain ----
  st_lv2_gain <-
    eq3 %>%
      dplyr::filter(Gtj > St) %>%
      dplyr::group_by(To) %>%
      dplyr::summarise(
        Gain = dplyr::n(),
        N = length(unique(eq3$Period)),
        Stationarity = "Active Gain",
        Test = ifelse(Gain == N, "Y", "N")
      ) %>%
      rbind(
        eq3 %>%
          dplyr::filter(Gtj < St) %>%
          dplyr::group_by(To) %>%
          dplyr::summarise(
            Gain = dplyr::n(),
            N = length(unique(eq3$Period)),
            Stationarity = "Dormant Gain",
            Test = ifelse(Gain == N, "Y", "N")
          )
      )

  # ---- Categorical loss ----
  st_lv2_loss <-
    eq4 %>%
      dplyr::filter(Lti > St) %>%
      dplyr::group_by(From) %>%
      dplyr::summarise(
        Loss = dplyr::n(),
        N = length(unique(eq4$Period)),
        Stationarity = "Active Loss",
        Test = ifelse(Loss == N, "Y", "N")
      ) %>%
      rbind(
        eq4 %>%
          dplyr::filter(Lti < St) %>%
          dplyr::group_by(From) %>%
          dplyr::summarise(
            Loss = dplyr::n(),
            N = length(unique(eq4$Period)),
            Stationarity = "Dormant Loss",
            Test = ifelse(Loss == N, "Y", "N")
          )
      )

  # Combine outputs into list
  intensity_tables <- list(lulc_table = lulc,
                           interval_lvl = level01,
                           category_lvlGain = list(categoryData = eq3,
                                                   categoryStationarity = st_lv2_gain),
                           category_lvlLoss = list(categoryData = eq4,
                                                   categoryStationarity = st_lv2_loss)
  )
  return(intensity_tables)
}

# close function

#' simulation_intensity_analysis: 
#' 
#' Perform intensity analysis (IA) on the time series of simulated LULC maps
#'
#' @param folder Folder with simulation maps
#' @param save_folder Folder to save IA results
#' @param base_name Base string for naming result file
#' @param scenario_name Name of scenario
#' @param lulc_agg LULC ratser attribute table
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

#dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

simulation_intensity_analysis <- function(
  folder, save_folder, base_name, scenario_name, lulc_agg
) {

  # Empty vector to store warnings
  warnings <- c()

  # Get all maps - tif or rds
  map_paths <- list.files(folder, full.names = TRUE)
  map_paths <- map_paths[grepl(".tif|.rds|.grd", map_paths)]
  # filter for 'simulated_LULC_simID_{simid}_year_2030.tif'
  map_paths <- map_paths[grepl("simulated_LULC_", map_paths)]
  # sort after years
  map_paths <- map_paths[order(as.numeric(stringr::str_match(basename(map_paths), "_(\\d{4})(_|\\.)")[, 2]))]
  print(map_paths)
  # If no maps found, return
  if (length(map_paths) == 0) {
    cat(paste0("No maps found in ", folder, ". Skipping folder.\n"))
    return()
  }
  # If one map found then message that IA is not possible
  if (length(map_paths) == 1) {
    cat(paste0("Only one map found in ", folder, " so intensity analysis is not possible. Skipping folder.\n"))
    return()
  }

  # Match year in map name by _(\d{4})(_|\.) and extract first capture group
  map_years <- stringr::str_match(basename(map_paths), "_(\\d{4})(_|\\.)")[, 2]
  # Load maps as raster stack
  map_stack <- stack(map_paths)
  # Rename layers to _year convention required by OpenLand::contingencyTable()
  names(map_stack) <- paste0("_", map_years)

  # Test if any rasters are empty
  if (any(is.na(minValue(map_stack)))) {
    warnings <- c(warnings, paste0("Folder contains empty raster(s)"))
  }

  #calculate pixel-wise cumulative changes
  pixel_cum_changes <- acc_changes(map_stack)

  #check that the number of pixels changing 5 times or more is not more than
  if (any(pixel_cum_changes[pixel_cum_changes$Cumulative_changes >= 5,]$Npixels > 1000)) {
    warnings <- c(warnings, paste0("The simulation contains more than 1000 pixels that are changing 5 times or more over the time steps"))
  }

  # Load LULC aggregation table
  lulc_agg <- readxl::read_xlsx(lulc_agg)
  #create a raster attribute table
  lulc_rat <- data.frame(ID = sort(unique(values(map_stack[[1]]))))
  # lulc_rat$lulc_name <- unlist(sapply(lulc_rat$ID, function(y) unique(unlist(lulc_agg[lulc_agg$Aggregated_ID == y, "Class_abbreviation"])), simplify = TRUE))

  # Perform intensity analysis
  Sim_IA <- intensityAnalysis(
    input_rasters = map_stack,
    pixelresolution = 100
  )

  # Add From/RTo class names to SIM_IA$lulc_table
  # Sim_IA$lulc_table$From_name <- unlist(sapply(Sim_IA$lulc_table$From, function(y) unique(unlist(lulc_rat[lulc_rat$ID == y, "lulc_name"])), simplify = TRUE))
  # Sim_IA$lulc_table$To_name <- unlist(sapply(Sim_IA$lulc_table$To, function(y) unique(unlist(lulc_rat[lulc_rat$ID == y, "lulc_name"])), simplify = TRUE))

  # Check that all values in SIM_IA$interval_lvl$PercentChange are greater than 0 and less than 5
  if (any(Sim_IA$interval_lvl$PercentChange == 0)) {
    warnings <- c(warnings, "Warning: Some values of total areal LULC class change between the time steps are 0%, Check the data.")
  }
  if (any(Sim_IA$interval_lvl$PercentChange > 5)) {
    warnings <- c(warnings, "Warning: Some values of total areal  LULC class change between the time steps are greater than 5%, Check the data.")
  }

  # Append pixel-wise cumulative changes to SIM_IA
  Sim_IA$pixel_cum_changes <- pixel_cum_changes

  # Append scenario name to base name
  base_name <- paste0(base_name, "_", scenario_name)
  # Assure folder exists
  dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

  # Save IA results as rds
  saveRDS(Sim_IA, file.path(save_folder, paste0(base_name, ".rds")))
  cat("Saved Intensity Analysis results as ", file.path(save_folder, paste0(base_name, ".rds")), "\n")

  return(warnings)
}

#' Perform Intensity Analysis (IA) for all simulations
#'
#' Loops over all subfolders in folder and calls simulation_intensity_analysis()
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
#' @param parallel Whether to run in parallel
#' @param lulc_agg LULC raster attribute table
#'
#' @docType methods

check_lulc_simulations <- function(
  folder, save_folder, base_name, parallel = FALSE, lulc_agg
) {
  if (length(list.dirs(config$InputDir, recursive = FALSE)) == 0) {
    # Process only this folder
    cat(paste0("Performing Intensity Analysis only for ", folder, " ...\n"))
    simulation_intensity_analysis(
      folder = folder,
      save_folder = save_folder,
      base_name = base_name,
      scenario_name = basename(folder)
    )
    return()
  }
  # Get all subfolders
  subfolders <- list.dirs(folder, recursive = FALSE, full.names = TRUE)
  # For each subfolder
  if (parallel == FALSE) {
    cat(paste0("Performing Intensity Analysis for ", length(subfolders),
               " subfolders in ", folder, " ...\n"))
    progress <- utils::txtProgressBar(
      min = 0,
      max = length(subfolders),
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
    cat(paste0("Performing Intensity Analysis for ", length(subfolders),
               " subfolders in ", folder, " in parallel with ", n_workers,
               " workers...\n"))
    future::plan(future::multisession, workers = n_workers)
  }
  # Shuffle `map_paths` to distribute memory load
  Simulation_warnings <- future.apply::future_sapply(
    # nolint start: indentation_linter
    subfolders,
    function(subfolder) {
      cat(paste0("Performing Intensity Analysis for ", subfolder, "...\n"))
      # Perform intensity analysis
      sim_warnings <- simulation_intensity_analysis(
        folder = subfolder,
        save_folder = save_folder,
        base_name = base_name,
        scenario_name = basename(subfolder),
        lulc_agg = lulc_agg
      )
      if (parallel == FALSE) {
        # Update progress bar
        utils::setTxtProgressBar(
          progress,
          value = utils::getTxtProgressBar(progress) + 1
        )
      }
      return(sim_warnings)
    })

  # Print message if warnings contain anything other than 'No warnings'
  if (any(Simulation_warnings != "No warnings.")) {
    cat("Warnings found in intensity analysis of LULC simulations, printing warnings:\n")
    print(Simulation_warnings)
    # Save warnings
    saveRDS(Simulation_warnings, file.path(config$OutputDir, "Simulation_warnings.rds"))
  } else {
    cat("No warnings found in intensity analysis of LULC simulations, results have been saved anyway\n")
  }
  
  # nolint end: indentation_linter
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
config <- config$CheckLULCC # only FocalLULCC

# if InputDir is not set, use bash variable
# FUTURE_EI_OUTPUT_DIR/LULCC_CH_OUTPUT_BASE_DIR
if (is.null(config$InputDir) || config$InputDir == "") {
  config$InputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR, bash_vars$LULCC_CH_OUTPUT_BASE_DIR)
}
# if OutputDir is not set, use bash variable CHECK_LULCC_OUTPUT_DIR
if (is.null(config$OutputDir) || config$OutputDir == "") {
  config$OutputDir <- file.path(bash_vars$FUTURE_EI_OUTPUT_DIR, bash_vars$CHECK_LULCC_OUTPUT_DIR)
}

# Check if InputDir is now set
if (is.null(config$InputDir) || config$InputDir == "") {
  stop("InputDir nor LULCC_CH_OUTPUT_BASE_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if InputDir exists
if (!dir.exists(config$InputDir)) {
  stop(cat("InputDir does not exist: ", config$InputDir, "\n"))
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
  stop("OutputDir nor CHECK_LULCC_OUTPUT_DIR set in FUTURE_EI_CONFIG_FILE")
}
# Check if OutputDir exists and create if not
if (!dir.exists(config$OutputDir)) {
  dir.create(config$OutputDir, recursive = TRUE)
}
# Check if BaseName is set
if (is.null(config$BaseName) || config$BaseName == "") {
  stop("BaseName not set in FUTURE_EI_CONFIG_FILE")
}
# Check if LULCC_M_CLASS_AGG is set
if (is.null(bash_vars$LULCC_M_CLASS_AGG) || bash_vars$LULCC_M_CLASS_AGG == "") {
  stop("LULCC_M_CLASS_AGG not set in FUTURE_EI_CONFIG_FILE")
}

cat("Performing Intensity Analysis on LULC simulations\n")

cat("Working directory is:", getwd(), "\n")
cat("Input directory set to:", config$InputDir, "\n")
cat("Output directory set to:", config$OutputDir, "\n")
cat("Base name set to:", config$BaseName, "\n")

t0 <- Sys.time()

# Perform intensity analysis on all simulations
check_lulc_simulations(
  folder = config$InputDir,
  save_folder = config$OutputDir,
  base_name = config$BaseName,
  parallel = config$Parallel,
  lulc_agg = file.path(bash_vars$LULCC_CH_HPC_DIR, bash_vars$LULCC_M_CLASS_AGG)
)

cat("Done performing itnensity analysis for LULC simulations!\n")
t0 <- Sys.time() - t0
cat("Time elapsed: ", t0, attr(t0, "units"), "\n")
