#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"

# Script to download, build and distribute the Docker image for the LULCC step
# Supports Docker, and Apptainer (Singularity) as a fallback only for downloading

# Assure LULCC_DOCKER_NAMESPACE, LULCC_DOCKER_REPO, and LULCC_DOCKER_VERSION are set
if [ -z "$LULCC_DOCKER_NAMESPACE" ] || [ -z "$LULCC_DOCKER_REPO" ] || [ -z "$LULCC_DOCKER_VERSION" ]; then
    log error "Please set the Docker Hub namespace, repository, and version of the LULCC Docker image. You can do this by setting the variables LULCC_DOCKER_NAMESPACE, LULCC_DOCKER_REPO, and LULCC_DOCKER_VERSION in the src/config.yml file."
    return
fi
# Docker Hub repository
namespace="$LULCC_DOCKER_NAMESPACE"  # can later be changed to an organizational account
repo="$LULCC_DOCKER_REPO"
version="$LULCC_DOCKER_VERSION"

# Check if the Docker or apptainer command is available
if command -v docker &> /dev/null; then
    # List docker images of the local machine for the repository, if any
    if docker image ls | grep "$repo" &> /dev/null; then
        log info "Available Docker images for $repo:"
        docker image ls | grep "$repo"
    fi

    # Check if the Docker image of the repository and version exists, else ask for download or build
    if docker image inspect "$repo:$version" &> /dev/null; then # default no
        log info "Docker image $repo:$version already exists. Want to rebuild it? [y/N]"
        read -r answer
        if [ "$answer" == "y" ]; then
            answer="b"
        else
            answer="s"
        fi
    else
        log info "Docker image $repo:$version does not exist. You can download it from Docker Hub or build it yourself. Do you want to build it? [(d)ownload/(b)uild/(S)kip]"
        read -r answer
    fi
    if [ "$answer" == "d" ]; then
        log info "Downloading Docker image $namespace/$repo:$version"
        docker pull "$namespace/$repo:$version"
    elif [ "$answer" == "b" ]; then
        log info "Building Docker image $repo:$version"
        docker build -t "$repo:$version" .
    elif [ "$answer" == "s" ]; then
        log info "Skipping build of Docker image $repo:$version"
    else
        log error "Invalid answer. Please choose (d)ownload/(b)uild/(s)kip"
        return
    fi

    # Check if the Docker image exists, else exit
    if ! docker image inspect "$repo:$version" &> /dev/null; then
        log error "Docker image $repo:$version does not exist. Please download or build it first."
        return
    fi

    # Check if the Docker image is available in the Docker Hub, else ask for push
    if docker image inspect "$namespace/$repo:$version" &> /dev/null; then
        log info "Docker image $namespace/$repo:$version already exists in Docker Hub. Skipping push."
    else
        log info "Docker image $namespace/$repo:$version does not exist in Docker Hub. Do you want to push it? [y/N]"
        read -r answer
        if [ "$answer" == "y" ]; then
            log info "Pushing Docker image $namespace/$repo:$version"
            docker tag "$repo:$version" "$namespace/$repo:$version"
            docker push "$namespace/$repo:$version"
        else
            log info "Skipping push of Docker image $namespace/$repo:$version"
        fi
    fi
elif command -v apptainer /dev/null; then
    log info "Apptainer is available from $(command -v apptainer)"
    # Create the Apptainer container directory if it does not exist
    if [ ! -d "$APPTAINER_CONTAINERDIR" ]; then
        log info "Creating Apptainer container directory: $APPTAINER_CONTAINERDIR"
        mkdir -p "$APPTAINER_CONTAINERDIR"
    else
        log info "Using existing Apptainer container directory: $APPTAINER_CONTAINERDIR"
    fi
    log info "Building SIF container $APPTAINER_CONTAINERDIR/${repo}_${version}.sif from docker base image $namespace/$repo:$version.
    Using submodules from $LULCC_CH_HPC_DIR/Model/Dinamica_models/LULCC_CH_ego_Submodels/*.ego"
    apptainer build \
      --build-arg "namespace=$namespace" --build-arg "repo=$repo" --build-arg "version=$version" \
      "$APPTAINER_CONTAINERDIR/${repo}_${version}.sif" "$SCRIPT_DIR/lulcc.def"
else
    log error "Neither Docker nor Apptainer is available. Please install one of them and make sure it is available in the PATH."
    return
fi

# To remove a specific image from the Docker Hub
# docker rmi $namespace/$repo:$version

# To remove all images of the repository from the local machine
# docker image rm $(docker image ls | grep $repo | awk '{print $3}')
