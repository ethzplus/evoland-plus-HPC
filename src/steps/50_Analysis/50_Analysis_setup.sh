#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Check LULCC Preparation

# Create conda environment
log info "Creating conda environment: ncp_analysis"
log debug "Using conda from: $CONDA_BIN"
log debug "Using requirement file: $SCRIPT_DIR/50_Analysis_env.yml"
$CONDA_BIN env create -f "$SCRIPT_DIR"/50_Analysis_env.yml

# Activate the conda environment
log debug "Activating conda environment: ncp_analysis"
source "$SCRIPT_DIR/../../de_activate.sh" ncp_analysis 1

# Not tracked in the environment.yml file:
# SSIMmap       : https://cloud.r-project.org/web/packages/SSIMmap/
log info "Installing SSIMmap R package"
R -e "pak::pkg_install(c('SSIMmap'), upgrade = FALSE)"
