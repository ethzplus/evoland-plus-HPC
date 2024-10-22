#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Check LULCC Preparation

# Create conda environment
log info "Creating conda environment: summarisation"
log debug "Using conda from: $CONDA_BIN"
log debug "Using requirement file: $SCRIPT_DIR/50_Summarisation_env.yml"
$CONDA_BIN env create -f "$SCRIPT_DIR"/50_Summarisation_env.yml

# Activate the conda environment
log debug "Activating conda environment: Summarisation"
source "$SCRIPT_DIR/../../de_activate.sh" Summarisation 1