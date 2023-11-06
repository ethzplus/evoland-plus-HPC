#!/bin/bash
# File for loading bash variables

# Try to call jq, if not load module
if ! (command -v jq &> /dev/null); then
    module load jq
    echo "Loaded jq module ($(jq --version))" || log error "Could not load jq module"
fi

# load bash_variables from config.json
config_file="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/config.json"
if [ -f "$config_file" ]; then
    eval "$(jq -r '.bash_variables | to_entries | .[] | "export \(.key)=\(.value)"' "$config_file")"
    echo "Loaded bash_variables from $config_file. Logging level: $FUTURE_EI_LOG_LEVEL"
else
    echo "No config.json file found"
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
fi

# If $TMPDIR not set, set it to $ALTERNATIVE_TMPDIR. Log it.
if [ -z "$TMPDIR" ]; then
    export TMPDIR="$ALTERNATIVE_TMPDIR"
    log debug "Setting TMPDIR to $TMPDIR"
else
    log debug "TMPDIR already set to $TMPDIR"
fi
