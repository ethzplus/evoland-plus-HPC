#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Focal Window Preparation

# Create conda environment
log info "Creating conda environment: focal_lulc"
log debug "Using conda from: $CONDA_BIN"
log debug "Using requirement file: $SCRIPT_DIR/20_FocalPrep_requirements.txt"
$CONDA_BIN env create -n focal_lulc -c conda-forge --file "$SCRIPT_DIR"/20_FocalPrep_requirements.txt

# Activate the conda environment
log debug "Activating conda environment: focal_lulc"
source "$SCRIPT_DIR/../../de_activate.sh" focal_lulc 1

# Export the conda environment
log debug "Exporting conda environment: focal_lulc"
$CONDA_BIN env export | grep -v "^prefix: " > "$SCRIPT_DIR"/20_FocalPrep_env.yml