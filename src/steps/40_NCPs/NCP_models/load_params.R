# Script to load the NCP parameters
#
# Load the parameters from environment variable NCP_PARAMS_YML if set,
# otherwise from ./40_NCPs_params.yml

# Load libraries
library(yaml)

# Try to load the parameters from the environment - catch error if not set
tryCatch(
  {
    # Load the parameters from env var NCP_PARAMS_YML if set, otherwise from
    # ./40_NCPs_params.yml
    if (Sys.getenv("NCP_PARAMS_YML") == "") {
      print(paste(
        "NCP_PARAMS_YML environment variable is not set. ",
        "Loading from ./40_NCPs_params.yml"
      ))
      params <- yaml.load_file(
        # Find config file in the same directory as the current script
        file.path(dirname(sys.frame(1)$ofile), "40_NCPs_params.yml")
      )
    } else {
      print(paste0(
        "Loading NCP parameters from ", Sys.getenv("NCP_PARAMS_YML")
      ))
      params <- yaml.load_file(Sys.getenv("NCP_PARAMS_YML"))
    }
  },
  error = function(err) {
    # If err is not a file not found error, re-throw the error
    if (!grepl("No such file or directory", err$message)) stop(err)
    # Throw error
    stop(
      paste0(
        "Loading the NCP parameters failed. Please set the NCP_PARAMS_YML ",
        "environment variable to the path of the 40_NCPs_params.yml file. ",
        "Got error: ", err$message
      )
    )
  }
)
