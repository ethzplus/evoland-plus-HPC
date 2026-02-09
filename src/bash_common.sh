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
    echo "WARNING: yq not found. Installation instructions: https://github.com/mikefarah/yq/#install"
    echo "Continuing without loading config.yml variables..."
    SKIP_YQ=1
fi

# load bash_variables from config.yml
config_file="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/config.yml"
if [ -f "$config_file" ] && [ -z "$SKIP_YQ" ]; then
    eval "$(yq e '.bash_variables | to_entries | .[] | "export " + .key + "=" + .value' "$config_file")"
    echo "Loaded bash_variables from $config_file. Logging level: $FUTURE_EI_LOG_LEVEL"
else
    if [ -n "$SKIP_YQ" ]; then
        echo "Skipping config.yml loading (yq not available)"
    else
        echo "No config.yml file found"
    fi
fi

# Add logging
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/logging.sh"
log debug "Loaded logging"

# Find out if micromamba > mamba > conda is installed and set bin path
export CONDA_BIN

# Check if MAMBA_EXE is set but points to non-existent file
if [ -n "$MAMBA_EXE" ] && [ ! -f "$MAMBA_EXE" ]; then
    log debug "MAMBA_EXE is set to '$MAMBA_EXE' but file does not exist. Will search for micromamba."
    unset MAMBA_EXE
fi

# If MAMBA_EXE not set or invalid, check multiple possible locations
if [ -z "$MAMBA_EXE" ]; then
    POSSIBLE_LOCATIONS=(
        "$HOME/.local/bin/micromamba"
        "/cluster/home/bblack/.local/bin/micromamba"
    )
    
    # Add custom location if set
    if [ -n "$MAMBA_EXE_CUSTOM" ]; then
        POSSIBLE_LOCATIONS+=("$MAMBA_EXE_CUSTOM")
    fi
    
    log debug "Searching for micromamba in: ${POSSIBLE_LOCATIONS[*]}"
    
    for location in "${POSSIBLE_LOCATIONS[@]}"; do
        log debug "Checking: $location"
        if [ -f "$location" ] && [ -x "$location" ]; then
            export MAMBA_EXE="$location"
            log debug "Found micromamba at: $MAMBA_EXE"
            break
        fi
    done
fi

if [ -n "$MAMBA_EXE" ] && [ -f "$MAMBA_EXE" ]; then
    # $MAMBA_EXE points to either micromamba or mamba
    CONDA_BIN=$MAMBA_EXE
    eval "$($CONDA_BIN shell hook -s bash)"
    log debug "Using mamba at $CONDA_BIN"
elif [ -x "$(command -v conda)" ]; then
    CONDA_BIN=$(which conda)
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
elif [ -f "$HOME/conda/bin/conda" ]; then
    # Check for conda in home directory
    CONDA_BIN="$HOME/conda/bin/conda"
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    # Check for miniconda3 in home directory
    CONDA_BIN="$HOME/miniconda3/bin/conda"
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
elif [ -f "$HOME/anaconda3/bin/conda" ]; then
    # Check for anaconda3 in home directory
    CONDA_BIN="$HOME/anaconda3/bin/conda"
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
else
    log error "Conda not found. Please ensure micromamba or conda is installed."
    log error "Checked micromamba locations: ${POSSIBLE_LOCATIONS[*]}"
    log error "Checked conda locations: \$HOME/conda/bin/conda, \$HOME/miniconda3/bin/conda, \$HOME/anaconda3/bin/conda"
    log error "Current MAMBA_EXE value: ${MAMBA_EXE:-not set}"
    log error "You can set MAMBA_EXE_CUSTOM environment variable to specify a custom location."
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
