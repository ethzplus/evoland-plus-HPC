#!/bin/bash
# The LULCC requires Dinamica EGO. Docker file will download the AppImage and add the
# path to the APP_DIR and DINAMICA_EGO_CLI variable in the Dockerfile.
# The AppImage can be downloaded from:
# https://dinamicaego.com/dinamica-7/
# or directly from:
# wget https://dinamicaego.com/nui_download/1792/

# Check if the APP_DIR and DINAMICA_EGO_CLI variables are set and exist
if [ -z "$APP_DIR" ] || [ -z "$DINAMICA_EGO_CLI" ]; then
    echo "APP_DIR or DINAMICA_EGO_CLI not set. Please set the variables in the .bash_variables file in the root directory of the repository."
    return
fi
if [ ! -f "$DINAMICA_EGO_CLI" ]; then
    echo "DINAMICA_EGO_CLI does not exist. Please set the variable in the .bash_variables file in the root directory of the repository."
    return
fi

# The external communication with R is described in
# https://dinamicaego.com/dinamica/dokuwiki/doku.php?id=external_communication
# Verion 1.0.4 of the R package is used. It can be downloaded from:
# https://dinamicaego.com/dinamica/dokuwiki/lib/exe/fetch.php?media=dinamica_1.0.4.tar.gz
if [ ! -f "$APP_DIR/dinamica_1.0.4.tar.gz" ]; then
    echo "Downloading Dinamica EGO R package"
    curl https://dinamicaego.com/dinamica/dokuwiki/lib/exe/fetch.php?media=dinamica_1.0.4.tar.gz -o "$APP_DIR/dinamica_1.0.4.tar.gz"
fi

# Create the conda environment
micromamba env create -n land_use -c conda-forge --file "$APP_DIR"/10_land_use_requirements.txt --yes

# Activate the conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate land_use

# Not tracked in the environment.yml file:
# Dinamica EGO connector - Install the Dinamica R package
R -e "install.packages('$APP_DIR/dinamica_1.0.4.tar.gz', repos = NULL, type = 'source')"

# Export the conda environment
micromamba env export | grep -v "^prefix: " > "$APP_DIR"/10_land_use_env.yml
