#####
#
#This script produces an RP indicator as a normalized aggregate (sum) of
#three landscape characteristics maps:
# 1.Degree of naturalness (DN): Calculated by applying a naturalness score from
#NCP_models/REC/lutable_naturality.csv to each LULC class within the simulated land-use map.
#2.Natural protected areas (NP): a binary map of 0=outside of protected areas, 1=inside of protected areas
#3. Water components (W):Calculated by computing inverse relative distance to lake coasts, 
#getting the highest value at lake coast and a decreasing value for 2km. 
# The output is a single map of Recreation potential.
#
#
#####

# #testing
# params <- list("data" = list(distlakes_path = "NULL",
#                              lulc = "C:/Users/bblack/switchdrive/C.3_Modelling/LULCC_CH_HPC/Results/Dinamica_simulated_LULC/BAU/v6/simulated_LULC_scenario_BAU_simID_v6_year_2060.tif"),
#                "REC" = list(lutable_nat_path = "C:/Users/bblack/switchdrive/Private/git_laptop/Future-EI/src/steps/40_NCPs/NCP_models/REC/lutable_naturality.csv"),
#                "run_params" = list(NCP_RUN_SCENARIO_ID = "1",
#                                    NCP_RUN_YEAR = "2020",
#                                    LULCC_M_EI_LAYER_DIR = "C:/Users/bblack/switchdrive/C.3_Modelling/LULCC_CH_HPC/Data/EI_intervention_layers")
# )

# Load libraries
library(terra)
library(yaml)

# Load the parameters into env by sourcing the ../load_params.R script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_dir <- dirname(sub(
  file_arg_name, "", initial_options[grep(file_arg_name, initial_options)]
))
source(file.path(script_dir, "..", "load_params.R"))
# Check all required parameters are set
if (is.null(params$data$distlakes_path)) {
  stop("params$data$distlakes_path is not set")
}
if (is.null(params$REC$lutable_nat_path)) {
  stop("params$REC$lutable_nat_path is not set")
}
if (is.null(params$run_params$NCP_RUN_SCENARIO_ID)) {
  stop("params$run_params$NCP_RUN_SCENARIO_ID is not set")
}
if (is.null(params$run_params$NCP_RUN_YEAR)) {
  stop("params$run_params$NCP_RUN_YEAR is not set")
}
if (is.null(params$run_params$LULCC_M_EI_LAYER_DIR)) {
  stop("params$run_params$LULCC_M_EI_LAYER_DIR is not set")
}
if (is.null(params$data$lulc)) {
  stop("params$data$lulc is not set")
}

# Vector output directory and output file name
results <- file.path(
  params$run_params$NCP_RUN_OUTPUT_DIR,
  params$run_params$NCP_RUN_SCENARIO_ID, "REC"
)
dir.create(results, showWarnings = FALSE, recursive = TRUE)
name_out <- paste("REC_S_CH_", params$run_params$NCP_RUN_YEAR, ".tif", sep = "")


# 1) Degree of naturalness ---------------------------------------------------

# Load current simulation/time step LULC 
lulc <- rast(params$data$lulc)

# Load the LULC naturality lookup table
lutable_NAT <-read.csv(params$REC$lutable_nat_path, header=T, sep=",") 

# Seperate LULC and HABITAT columns and convert to matrix
m <- as.matrix(lutable_NAT[, c("LULC", "HABITAT")])

# Reclassify the LULC raster
nc<-classify(lulc,m1)


# 2) Natural/protected areas -------------------------------------------------

Current_PA_path <- list.files(file.path(params$run_params$LULCC_M_EI_LAYER_DIR, "Future_EI",
                                        paste0("EI_ID",params$run_params$NCP_RUN_SCENARIO_ID),
                                        params$run_params$NCP_RUN_YEAR),
                                        pattern = ".tif",
                                        full.names = TRUE)

# Exclude tif.ovr files
Current_PA_path <- Current_PA_path[!grepl(".ovr", Current_PA_path)]

# Load the protected areas shapefile
pa<-  terra::rast(Current_PA_path)

# Reclassify the raster such that NA values are 0
# This gives a raster of 0=outside of protected areas, 1=inside of protected areas
m<-c(NA,0)
m1<-matrix(m,byrow=T, ncol=2)
ppa <- terra::classify(pa,m1)

# 3.) Water component  -----------------------------------------------------

# Load raster with a 2km  distance buffer around lakes (done in arcgis Pro)
distlakes<- terra::rast(params$data$distlakes_path)  

# Processing distance layer to give higher value close to the lake
# Normalizing values
nx <- minmax(distlakes)    
rn <- (distlakes - nx[1,]) / (nx[2,] - nx[1,])

# Inverting values 
rn2<- rn
values(rn2)<-1-values(rn)

# Classifying NAs to 0
m<-c(NA, 0)
m1<-matrix(m, byrow=T, ncol=2)
rn2<-classify(rn2, m1)

# Mask layer to the extent of the lulc raster
rn2<-mask(rn2, lulc)

# Extend and crop the raster to the project extent
rn2<-terra::extend(rn2,ext(params$proj$ext))
rn2<-terra::crop(rn2,ext(params$proj$ext))


# 4) Merge components and export  --------------------------------------------------

# Combine the three components
nat <- ppa+nc+rn2

# Normalize the values to 0-1
nat2 <- nat/3 

# Export
writeRaster(nat2, file.path(results, name_out), overwrite = TRUE)
print(paste("REC layer created: ", file.path(results, name_out)))
