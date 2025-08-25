#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"

# Create Python environment
log info "Creating conda environment: ncps-python"
$CONDA_BIN env create -f "$SCRIPT_DIR"/40_NCPs_python_env.yml -y

# Create R environment  
log info "Creating conda environment: ncps-r"
$CONDA_BIN env create -f "$SCRIPT_DIR"/40_NCPs_r_env.yml -y

# Install R packages in R environment
log debug "Activating conda environment: ncps-r"
source "$SCRIPT_DIR/../../de_activate.sh" ncps-r 1
log info "Installing meteor and ecocrop R packages"
R -e "pak::pkg_install(c('meteor', 'cropmodels/Recocrop'), upgrade = FALSE)"

log info "Both environments created successfully"