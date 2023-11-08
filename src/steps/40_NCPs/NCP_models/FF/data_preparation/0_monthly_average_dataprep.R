###################
#
# The aim of this script is to prepare the monthly average data for the ecocrop model
#
# The script reads precipitation and temperature datasets and preprocesses them.
#
###################

library(stringr)
library(terra)
library(raster)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "..", "40_NCPs_params.yml")
)

# local

ecocrop_rawdata_dir <- params$FF$ecocrop_rawdata_dir

scratch <- paste(ecocrop_rawdata_dir, "scratch", sep = "/")

pr <- paste(ecocrop_rawdata_dir, "prec", sep = "/")
pavg_dir <- params$FF$pavg_dir

ta <- paste(ecocrop_rawdata_dir, "temp", sep = "/")
tavg_dir <- params$FF$tavg_dir

#----------------------------------1992-1997 - PRECIPITATION

#loading and renaming

list_names_13 <- list()

for (i in 1:length(list.files(pr))) {
  a <- list.files(pr)[i]
  load(paste(pr, a, sep = "/"))
  name <- substring(a, 1, nchar(a) - 6)
  assign(name, X)
  list_names_13 <- append(list_names_13, name)

  print(paste(round(i * 100 / length(list.files(pr)), 1), "%", sep = " "))
}

# monthly means 

list_months <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12")

for (i in 1:length(list_months)) {

  mnth <- list_months[i]
  a <- grep(paste(mnth, "\\b", sep = ""), list_names_13)

  b <- list_names_13[a]

  c <- mean(get(paste(b[1])), get(paste(b[2])), get(paste(b[3])), get(paste(b[4])), get(paste(b[5])))
  ex_n <- paste("PrecAvg_13-18_", mnth, ".tif", sep = "")

  c <- terra::rast(c)
  terra::writeRaster(c, paste(pavg_dir, ex_n, sep = "/"), overwrite = TRUE)

  print(paste(round(i * 100 / length(list_months), 1), "%", sep = " "))

}
#---------------------------------- - AVG Temp

#loading and renaming

list_names_13 <- list()

for (i in 1:length(list.files(ta))) {
  a <- list.files(ta)[i]
  load(paste(ta, a, sep = "/"))
  name <- substring(a, 1, nchar(a) - 6)
  assign(name, X)
  list_names_13 <- append(list_names_13, name)

  print(paste(round(i * 100 / length(list.files(ta)), 1), "%", sep = " "))
}

# monthly means 

list_months <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12")

for (i in 1:length(list_months)) {

  mnth <- list_months[i]
  a <- grep(paste(mnth, "\\b", sep = ""), list_names_13)

  b <- list_names_13[a]

  c <- mean(get(paste(b[1])), get(paste(b[2])), get(paste(b[3])), get(paste(b[4])), get(paste(b[5])))
  ex_n <- paste("TavgAvg_13-18_", mnth, ".tif", sep = "")

  c <- terra::rast(c)
  terra::writeRaster(c, paste(tavg_dir, ex_n, sep = "/"), overwrite = TRUE)

  print(paste(round(i * 100 / length(list_months), 1), "%", sep = " "))

}


rm(list = ls())

