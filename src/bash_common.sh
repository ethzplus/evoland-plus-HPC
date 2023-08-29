# File for loading bash variables

# load bash_variables from config.json
config_file="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/config.json"
if [ -f $config_file ]; then
    eval $(jq -r '.bash_variables | to_entries | .[] | "export \(.key)=\(.value)"' $config_file)
    echo "Loaded bash_variables from $config_file. Logging level: $FUTURE_EI_LOG_LEVEL"
else
    echo "No config.json file found"
fi

# Add logging
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/logging.sh"
log debug "Loaded logging"

# Find out if micromamba > miniconda > conda is installed and set bin path
if [ -f $MAMBA_EXE ]; then
    # points to either micromamba or mamba
    export CONDA_BIN=$MAMBA_EXE
    eval "$($CONDA_BIN shell hook -s bash)"
    log debug "Using mamba at $CONDA_BIN"
elif [ -x "$(command -v conda)" ]; then
    export CONDA_BIN=$(which conda)
    eval "$($CONDA_BIN shell.bash hook)"
    log debug "Using conda at $CONDA_BIN"
else
    log error "Conda not found"
fi
