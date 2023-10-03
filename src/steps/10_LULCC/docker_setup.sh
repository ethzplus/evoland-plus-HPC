#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../bash_common.sh"

# Script to download, build and distribute the Docker image for the LULCC step
# Supports Docker, and Apptainer (Singularity) as a fallback only for downloading

# Docker Hub repository
version="0.1"
namespace="cbueth"  # can later be changed to an organizational account
repo="lulcc"

apptainer=$(command -v singularity)  #$(command -v apptainer)

# Check if the Docker or apptainer command is available
if command -v docker &> /dev/null; then
    # List docker images of the local machine for the repository, if any
    if docker image ls | grep $repo &> /dev/null; then
        log info "Available Docker images for $repo:"
        docker image ls | grep $repo
    fi

    # Check if the Docker image of the repository and version exists, else ask for download or build
    if docker image inspect $repo:$version &> /dev/null; then
        log info "Docker image $repo:$version already exists. Skipping build."
    else
        log info "Docker image $repo:$version does not exist. You can download it
         from Docker Hub or build it yourself. Do you want to build it? [
         (d)ownload/(b)uild/(s)kip]"
        read -r answer
        if [ "$answer" == "d" ]; then
            log info "Downloading Docker image $namespace/$repo:$version"
            docker pull $namespace/$repo:$version
        elif [ "$answer" == "b" ]; then
            log info "Building Docker image $repo:$version"
            docker build -t $repo:$version .
        elif [ "$answer" == "s" ]; then
            log info "Skipping build of Docker image $repo:$version"
        else
            log error "Invalid answer. Please choose (d)ownload/(b)uild/(s)kip"
            return
        fi
    fi

    # Check if the Docker image exists, else exit
    if ! docker image inspect $repo:$version &> /dev/null; then
        log error "Docker image $repo:$version does not exist. Please download or build it first."
        return
    fi

    # Check if the Docker image is available in the Docker Hub, else ask for push
    if docker image inspect $namespace/$repo:$version &> /dev/null; then
        log info "Docker image $namespace/$repo:$version already exists in Docker Hub. Skipping push."
    else
        log info "Docker image $namespace/$repo:$version does not exist in Docker Hub. Do you want to push it? [y/n]"
        read -r answer
        if [ "$answer" == "y" ]; then
            log info "Pushing Docker image $namespace/$repo:$version"
            docker tag $repo:$version $namespace/$repo:$version
            docker push $namespace/$repo:$version
        else
            log info "Skipping push of Docker image $namespace/$repo:$version"
        fi
    fi
elif command -v "$apptainer" /dev/null; then
    log info "Apptainer is available from $(command -v singularity)"
    log info "Downloading Docker image $namespace/$repo:$version"
    $apptainer pull docker://$namespace/$repo:$version
else
    log error "Neither Docker nor apptainer is available. Please install one of them and make sure it is available in the PATH."
    return
fi

# To remove a specific image from the Docker Hub
# docker rmi $namespace/$repo:$version

# To remove all images of the repository from the local machine
# docker image rm $(docker image ls | grep $repo | awk '{print $3}')
