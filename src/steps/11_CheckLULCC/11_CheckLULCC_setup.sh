#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Check LULCC Preparation

# Create conda environment
log info "Creating conda environment: check_lulc"
log debug "Using conda from: $CONDA_BIN"
log debug "Using requirement file: $SCRIPT_DIR/11_checklulcc_env.yml"
$CONDA_BIN env create -f "$SCRIPT_DIR"/11_checklulcc_env.yml

# Activate the conda environment
log debug "Activating conda environment: check_lulc"
source "$SCRIPT_DIR/../../de_activate.sh" check_lulc 1