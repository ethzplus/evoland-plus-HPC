#!/bin/bash
# File for loading bash variables

set -e  # Stop script if any command fails in parent or child scripts
set -o pipefail  # Return exit status of the last command in the pipe that failed

# Make sure only sourced once
if [ -n "$BASH_COMMON_LOADED" ]; then
    return
else
    BASH_COMMON_LOADED=1
fi

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

preparation_script="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/steps/00_Preparation.sh"
# shellcheck source=src/steps/00_Preparation.sh
source "$preparation_script"
