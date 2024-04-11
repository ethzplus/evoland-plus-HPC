#!/bin/bash
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/bash_common.sh"
# Script to de-/activate environment passed as first argument
# Only de/activate depending on second argument: 1 - activate, 0 - deactivate

# if $1 not set do nothing
if [ -z "$1" ]; then
    echo "No environment passed as argument. Aborting."
    exit 1
fi

# if $CONDA_BIN not set exit
if [ -z "$CONDA_BIN" ]; then
    echo "CONDA_BIN not set. Aborting."
    exit 1
fi

# Is any environment already active and $2 is None or 0?
if [ -n "$CONDA_DEFAULT_ENV" ] && [ -z "$2" ] || [ "$2" == "0" ]; then
    # Deactivate environment
    if [[ $CONDA_BIN == *"micromamba"* ]]; then
        echo "Deactivating environment $CONDA_DEFAULT_ENV with micromamba"
        micromamba deactivate
    elif [[ $CONDA_BIN == *"mamba"* ]]; then
        echo "Deactivating environment $CONDA_DEFAULT_ENV with mamba"
        mamba deactivate
    elif [[ $CONDA_BIN == *"conda"* ]]; then
        echo "Deactivating environment $CONDA_DEFAULT_ENV with conda"
        conda deactivate
    else
        echo "No conda installation found"
        exit 1
    fi
# Only if $CONDA_DEFAULT_ENV is not set and $2 is None or 1
elif [ -z "$CONDA_DEFAULT_ENV" ] && [ -z "$2" ] || [ "$2" == "1" ]; then
    # Activate environment
    # Add bash hook - different for micromamba, mamba and conda
    # if $CONDA_BIN end with micromamba
    if [[ $CONDA_BIN == *"micromamba"* ]]; then
        echo "Activating environment $1 with micromamba"
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate "$1"
    elif [[ $CONDA_BIN == *"mamba"* ]]; then
        echo "Activating environment $1 with mamba"
        eval "$(mamba shell hook -s bash)"
        mamba activate "$1"
    elif [[ $CONDA_BIN == *"conda"* ]]; then
        echo "Activating environment $1 with conda"
        eval "$(conda shell.bash hook)"
        conda activate "$1"
    else
        echo "No conda installation found"
        exit 1
    fi
fi