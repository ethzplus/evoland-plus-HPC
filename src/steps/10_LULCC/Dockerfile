# Docker image LULCC uses
# https://github.com/cbueth/dinamica-ego-docker/
# Unpublished image - build locally with dinamica-ego-docker
# hadolint ignore=DL3006
FROM dinamica-ego
LABEL authors="Carlson Büth" \
      description="Docker image for the LULCC (land use land cover change) model"

ENV DINAMICA_EGO_CLI_LOG_LEVEL=4
# Change log level by passing --env DINAMICA_EGO_CLI_LOG_LEVEL=7 when running
# UNCONDITIONAL=0, ERROR=1, WARNING=2, RESULT=3,
# INFO=4, INFO2=5, DEBUG=6, DEBUG2=7

# Check if $DINAMICA_EGO_CLI and $APP_DIR are set and exist
RUN if [ -z "$DINAMICA_EGO_CLI" ]; then echo "DINAMICA_EGO_CLI not found"; exit 1; fi \
 && if [ -z "$APP_DIR" ]; then echo "APP_DIR not found"; exit 1; fi \
 && if [ ! -f "$DINAMICA_EGO_CLI" ]; then echo "DINAMICA_EGO_CLI does not exist"; exit 1; fi \
 && if [ ! -d "$APP_DIR" ]; then echo "APP_DIR does not exist"; exit 1; fi

# hadolint ignore=DL3009
RUN apt-get update && apt-get install libcurl4-openssl-dev=8.5.0-2 \
    libssl-dev=3.1.5-1 libfreetype-dev=2.13.2+dfsg-1+b1 -y --no-install-recommends

RUN install2.r --error --skipinstalled \
    curl pak pkgdepends
# Use pak to install all R and system dependencies
COPY 10_land_use_requirements.txt $APP_DIR/10_land_use_requirements.txt
RUN R -e "pak::pkg_install(readLines(\"$APP_DIR/10_land_use_requirements.txt\"), upgrade = TRUE)" \
 && rm -rf /var/lib/apt/lists/*
# Show dependency tree with pak::pkg_deps_tree(readLines("10_land_use_requirements.txt"))

WORKDIR $MODEL_DIR

# Define entrypoint for the DinamicaConsole in the AppImage, within Micromamba environment
ENTRYPOINT ["/bin/bash", "-c", "$DINAMICA_EGO_CLI -disable-parallel-steps -log-level $DINAMICA_EGO_CLI_LOG_LEVEL $MODEL_DIR/Dinamica_models/LULCC_CH.ego"]

