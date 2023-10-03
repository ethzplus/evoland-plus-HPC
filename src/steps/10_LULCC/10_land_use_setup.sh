#!/bin/bash
# ---> Rather use containerized Docker image,
#      see src/steps/10_LULCC/docker_setup.sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"
# Setup script for the land use conda environment

# The LULCC requires Dinamica EGO. Download the AppImage and add the path to the
# DINAMICA_EGO_DIR and DINAMICA_EGO_APP variable in the .bash_variables file in the
# root directory of the repository. The AppImage can be downloaded from:
# https://dinamicaego.com/dinamica-7/
# or directly from:
# wget https://dinamicaego.com/nui_download/1792/

# Check if the DINAMICA_EGO_DIR and DINAMICA_EGO_APP variables are set and exist
if [ -z "$DINAMICA_EGO_DIR" ] || [ -z "$DINAMICA_EGO_APP" ]; then
    log error "DINAMICA_EGO_DIR or DINAMICA_EGO_APP not set. Please set the variables in the .bash_variables file in the root directory of the repository."
    return
fi
if [ ! -f "$DINAMICA_EGO_APP" ]; then
    log error "DINAMICA_EGO_APP does not exist. Please set the variable in the .bash_variables file in the root directory of the repository."
    return
fi

# The external communication with R is described in
# https://dinamicaego.com/dinamica/dokuwiki/doku.php?id=external_communication
# Verion 1.0.4 of the R package is used. It can be downloaded from:
# https://dinamicaego.com/dinamica/dokuwiki/lib/exe/fetch.php?media=dinamica_1.0.4.tar.gz
if [ ! -f "$DINAMICA_EGO_DIR/dinamica_1.0.4.tar.gz" ]; then
    log debug "Downloading Dinamica EGO R package"
    wget https://dinamicaego.com/dinamica/dokuwiki/lib/exe/fetch.php?media=dinamica_1.0.4.tar.gz -O "$DINAMICA_EGO_DIR/dinamica_1.0.4.tar.gz"
fi

# Create the conda environment
log info "Creating conda environment: land_use"
log debug "Using conda from: $CONDA_BIN"
log debug "Using requirement file: $SCRIPT_DIR/10_land_use_requirements.txt"
$CONDA_BIN create -n land_use -c conda-forge --file "$SCRIPT_DIR"/10_land_use_requirements.txt

# Activate the conda environment
log debug "Activating conda environment: land_use"
source "$SCRIPT_DIR/../../de_activate.sh" land_use 1

# Not tracked in the environment.yml file:
# Dinamica EGO connector - Install the Dinamica R package
log info "Installing Dinamica EGO R package"
R -e "install.packages('$DINAMICA_EGO_DIR/dinamica_1.0.4.tar.gz', repos = NULL, type = 'source')"

# Export the conda environment
log debug "Exporting conda environment: land_use"
$CONDA_BIN env export | grep -v "^prefix: " > "$SCRIPT_DIR"/10_land_use_env.yml
