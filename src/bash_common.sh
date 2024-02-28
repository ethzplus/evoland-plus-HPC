#!/bin/bash
# File for loading bash variables

# Try to call yq, if not refer to installation instructions (https://github.com/mikefarah/yq/#install)
if ! (command -v yq &> /dev/null); then
    echo "yq not found. Installation instructions: https://github.com/mikefarah/yq/#install"
fi

# load bash_variables from config.yml
config_file="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/config.yml"
if [ -f "$config_file" ]; then
    eval "$(yq e '.bash_variables | to_entries | .[] | "export " + .key + "=" + .value' "$config_file")"
    echo "Loaded bash_variables from $config_file. Logging level: $FUTURE_EI_LOG_LEVEL"
else
    echo "No config.yml file found"
fi

# Add logging
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/logging.sh"
log debug "Loaded logging"

# Find out if micromamba > mamba > conda is installed and set bin path
export CONDA_BIN
if [ -f "$MAMBA_EXE" ]; then
    # $MAMBA_EXE points to either micromamba or mamba
    CONDA_BIN=$MAMBA_EXE
    eval "$($CONDA_BIN shell hook -s bash)"
    log debug "Using mamba at $CONDA_BIN"
elif [ -x "$(command -v conda)" ]; then
    CONDA_BIN=$(which conda)
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
else
    log error "Conda not found"
    exit 1
fi

# If $TMPDIR not set, set it to $ALTERNATIVE_TMPDIR. Log it.
if [ -z "$TMPDIR" ]; then
    export TMPDIR="$ALTERNATIVE_TMPDIR"
    log debug "Setting TMPDIR to $TMPDIR"
else
    log debug "TMPDIR already set to $TMPDIR"
fi

# Use scratch for output folders
log debug "Using scratch for output folders:"
LULCC_CH_OUTPUT_BASE_DIR=$TMPDIR$LULCC_CH_OUTPUT_BASE_DIR
log debug "LULCC_CH_OUTPUT_BASE_DIR: $LULCC_CH_OUTPUT_BASE_DIR"
NCP_OUTPUT_BASE_DIR=$TMPDIR$NCP_OUTPUT_BASE_DIR
log debug "NCP_OUTPUT_BASE_DIR: $NCP_OUTPUT_BASE_DIR"
NSDM_OUTPUT_BASE_DIR=$TMPDIR$NSDM_OUTPUT_BASE_DIR
log debug "NSDM_OUTPUT_BASE_DIR: $NSDM_OUTPUT_BASE_DIR"
FOCAL_OUTPUT_BASE_DIR=$TMPDIR$FOCAL_OUTPUT_BASE_DIR
log debug "FOCAL_OUTPUT_BASE_DIR: $FOCAL_OUTPUT_BASE_DIR"
SPAGG_OUTPUT_BASE_DIR=$TMPDIR$SPAGG_OUTPUT_BASE_DIR
log debug "SPAGG_OUTPUT_BASE_DIR: $SPAGG_OUTPUT_BASE_DIR"
