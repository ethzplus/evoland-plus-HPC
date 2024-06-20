#!/bin/bash
# Prepare NCPs before running for every LULC map/year.
# The data this script prepares is used in the `run_all_ncps.sh` script.
# There are no passed arguments to this script, as the preparation is
# independent of the LULC map/year.
# Usage: bash prepare_ncps.sh
# Example: bash prepare_ncps.sh
# ------------------------------------------------------------------------------

# Script Variables
# ----------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../../bash_common.sh"
# Activate the conda environment with R and python packages
log debug "Activating conda environment: ncps"
source "$SCRIPT_DIR/../../../de_activate.sh" ncps 1
# Check if the conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "ncps" ]; then
    log error "The conda environment is not activated: ncps"
    exit 1
fi

source "$SCRIPT_DIR"/script_utils.sh

# ------------------------------------------------------------------------------

# Check if the required arguments are provided
if [ "$#" -ne 0 ]; then
    log error "Illegal number of parameters. Expected 0, Received $#"
    log error "Usage: bash prepare_ncps.sh"
    exit 1
fi

# ------------------------------------------------------------------------------

# Prepare NCPs

log info "Preparing NCPs"

## FF - Food and feed
## Indicator: Crop production potential (`ecocrop`)
log info "Preparing FF - Food and feed"
run_scripts "$SCRIPT_DIR/FF/data_preparation/0_FF_ecocrop.R"