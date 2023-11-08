#####
#
# The aim of this script is to reunite the different InVEST models for each
# region previously calculated in the other code files 1_CAR... & 2_CAR_...
# into one final map
#
#####

library(terra)
library(yaml)

# Load the parameters from ../40_NCPs_params.yml
params <- yaml.load_file(
  # Find config file relative to the location of the current script
  file.path(dirname(sys.frame(1)$ofile), "..", "40_NCPs_params.yml")
)

# Folder paths
setwd(params$CAR$wd)

dir18 <- params$CAR$dir18 # folder containing invest models
res18 <- params$CAR$res18 # folder to save results
n18 <- params$CAR$n18 # name of the final raster

# Function to apply for each time period

carbon_process_3 <- function(ind_dir, results, name_out) {

  #--- move all final outputs to new folder (this part scraps outputs from invest folders and thus copies the outputs. /!\ generates a copy of the files)
  # Lists all files in the input folder.
  list_files <- list.files(ind_dir)
  # Creates a new directory named 'tot_c_united' within the input folder.
  dir.create(paste(ind_dir, "tot_c_united", sep = "/"))

  print("directory created")
  # In this loop the individual carbon stok rasters are copied from their respective directory
  # to a unified directory 'tot_c_united'
  for (i in 1:length(list_files)) {

    print(list_files[i])
    path1 <- paste(ind_dir, list_files[i], sep = "/")
    name <- paste("tot_c_cur_", list_files[i], ".tif", sep = "")
    path2 <- paste(path1, name, sep = "/")
    file.copy(path2, paste(ind_dir, "tot_c_united", sep = "/"))
  }

  print("files copied")

  #--- bind together each file

  # creating empty raster to fill
  # TODO: Change for different countries
  bind <- rast()
  crs(bind) <- params$proj$crs
  ext(bind) <- c(params$proj$ext)
  res(bind) <- params$proj$res


  # merging each individual carbon stock raster into the empty raster created
  # above

  list_lu <- list.files(paste(ind_dir, "tot_c_united", sep = "/"))

  for (i in 1:length(list_lu)) {
    name <- list_lu[i]
    # loading a raster file from the folder
    nr <- rast(paste(ind_dir, "tot_c_united", name, sep = "/"))
    crs(nr) <- params$proj$crs
    # binding using the max function -> if there is an overlap, the max value is preserved
    bind <- mosaic(bind, nr, fun = "max")
    print(paste(name, " added"))
  }

  crs(bind) <- params$proj$crs
  bind <- extend(bind, c(params$proj$ext))

  writeRaster(bind, (paste(results, name_out, sep = "/")), overwrite = TRUE)

  print("final raster written")
}

#--- Applying for each period

carbon_process_3(dir18, res18, n18)
