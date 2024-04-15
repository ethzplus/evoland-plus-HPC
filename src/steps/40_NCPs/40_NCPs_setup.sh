#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Setup script for the ncps conda environment

# The NCPs require python and R packages, both are fetched from conda-forge

# Create the conda environment
log info "Creating conda environment: ncps"
log debug "Using requirement file: $SCRIPT_DIR/40_NCPs_env.yml"
# If $CONDA_BIN is conda print warning that it's solver might have problems with the environment.yml file
if [[ $CONDA_BIN == *conda ]]; then
    log warning "If conda has problems with the 40_NCPs_env.yml file, use micromamba instead, or add the libmamba solver (https://conda.github.io/conda-libmamba-solver/getting-started/)"
fi
$CONDA_BIN create -n ncps -c conda-forge --file "$SCRIPT_DIR"/40_NCPs_env.yml -y

# Activate the conda environment
log debug "Activating conda environment: ncps"
source "$SCRIPT_DIR/../../de_activate.sh" ncps 1

# Not tracked in the environment.yml file:
# meteor        : https://cran.r-project.org/web/packages/meteor/
# ecocrop       : https://github.com/cropmodels/Recocrop
log info "Installing meteor and ecocrop R packages"
R -e "pak::pkg_install(c('meteor', 'cropmodels/Recocrop'), upgrade = FALSE)"
