#!/bin/bash
# Land use HPC step for the batch job src/future-ei-pipeline.sh
# If called separately for testing, source the bash_common.sh script first

# Assure Apptainer (/Singularity) or Docker is available
if ! command -v apptainer &> /dev/null; then
    if ! command -v docker &> /dev/null; then
        log error "Neither Apptainer nor Docker is available. Please install one of them and make sure it is available in the PATH."
        return
    fi
fi
# Assure LULCC_CH_HPC_DIR is set
if [ -z "$LULCC_CH_HPC_DIR" ]; then
    log error "LULCC_CH_HPC_DIR is not set. Please set it to the directory of the LULCC-CH-HPC repository."
    return
fi
# Assure LULCC_DOCKER_IMAGE is set
if [ -z "$LULCC_DOCKER_IMAGE" ]; then
    log error "LULCC_DOCKER_IMAGE is not set. Please set it to the Docker image to use."
    return
fi
# Lazy load LULCC_DOCKER_IMAGE if not available
if ! apptainer exec docker://"$LULCC_DOCKER_IMAGE" true; then
    log info "Lazy loading Docker image $LULCC_DOCKER_IMAGE"
    apptainer pull docker://"$LULCC_DOCKER_IMAGE"
fi

# Docker image is available, run the container, preferably with Apptainer
log info "Running Docker image $LULCC_DOCKER_IMAGE with $LULCC_CH_HPC_DIR mounted to /model"
if ! command -v apptainer &> /dev/null; then
    log debug "Using Apptainer from $(command -v apptainer)"
    apptainer exec --bind "$LULCC_CH_HPC_DIR":/model docker://"$LULCC_DOCKER_IMAGE"
else
    log debug "Using docker from $(command -v docker)"
    docker run -v "$LULCC_CH_HPC_DIR":/model -it "$LULCC_DOCKER_IMAGE"
fi