#!/bin/bash
# Land use HPC step for the batch job src/future-ei-pipeline.sh
# If called separately for testing, source the bash_common.sh script first

# log SLURM variables
log debug "SLURM_JOB_ID: $SLURM_JOB_ID"
log debug "SLURM_JOB_NAME: $SLURM_JOB_NAME"
log debug "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
log debug "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
log debug "SLURM_JOB_CPUS_PER_NODE: $SLURM_JOB_CPUS_PER_NODE"

# Transfer simulation input parameters - prepare

# Assure Apptainer (/Singularity) or Docker is available
if ! (command -v apptainer &> /dev/null || command -v docker &> /dev/null); then
    log error "Neither Apptainer nor Docker is available. Please install one of them and make sure it is available in the PATH."
    return
fi
# Assure LULCC_CH_HPC_DIR is set
if [ -z "$LULCC_CH_HPC_DIR" ]; then
    log error "LULCC_CH_HPC_DIR is not set. Please set it to the directory of the LULCC-CH-HPC repository."
    return
fi
# Assure LULCC_DOCKER_NAMESPACE, LULCC_DOCKER_REPO, and LULCC_DOCKER_VERSION are set
if [ -z "$LULCC_DOCKER_NAMESPACE" ] || [ -z "$LULCC_DOCKER_REPO" ] || [ -z "$LULCC_DOCKER_VERSION" ]; then
    log error "Please set the Docker Hub namespace, repository, and version of the LULCC Docker image. You can do this by setting the variables LULCC_DOCKER_NAMESPACE, LULCC_DOCKER_REPO, and LULCC_DOCKER_VERSION in the src/config.yml file."
    return
fi
lulcc_docker_image="$LULCC_DOCKER_NAMESPACE/$LULCC_DOCKER_REPO:$LULCC_DOCKER_VERSION"

# Run the container, preferably with Apptainer
log info "Running Docker image $lulcc_docker_image with $LULCC_CH_HPC_DIR mounted to /model"
if command -v apptainer &> /dev/null; then
    log debug "Using Apptainer from $(command -v apptainer) with container $APPTAINER_CONTAINERDIR/${LULCC_DOCKER_REPO}_${LULCC_DOCKER_VERSION}.sif"
    log info "apptainer run --bind \"$FUTURE_EI_OUTPUT_DIR:/tmp,$LULCC_CH_HPC_DIR\":/model \"$APPTAINER_CONTAINERDIR/${LULCC_DOCKER_REPO}_${LULCC_DOCKER_VERSION}.sif\""
    apptainer run --bind "$FUTURE_EI_OUTPUT_DIR:/tmp,$LULCC_CH_HPC_DIR":/model "$APPTAINER_CONTAINERDIR/${LULCC_DOCKER_REPO}_${LULCC_DOCKER_VERSION}.sif"
else
    log debug "Using docker from $(command -v docker)"
    log warning "Docker support is experimental and may not work as expected."
    docker run -v "$LULCC_CH_HPC_DIR":/model -it "$lulcc_docker_image"
fi

# Check for ERROR in control table csv ($LULCC_M_SIM_CONTROL_TABLE)
log debug "Loading simulation control file: $LULCC_M_SIM_CONTROL_TABLE"
# Check if ERROR in any row of the control table
if grep -q ERROR "$LULCC_M_SIM_CONTROL_TABLE"; then
  # Log the number of occurrences of ERROR
  log error "Found $(grep -c ERROR "$LULCC_M_SIM_CONTROL_TABLE") occurrences of ERROR in simulation control file $LULCC_M_SIM_CONTROL_TABLE"
  return
else
  log info "No ERROR found in simulation control file $LULCC_M_SIM_CONTROL_TABLE"
fi
