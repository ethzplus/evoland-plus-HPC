# Docker image the LULCC use

# Use Ubuntu 22.10
FROM ubuntu:22.10
LABEL authors="Carlson Büth" \
      version="0.1" \
      description="Docker image for the LULCC (land use land cover change) model"

# Environment variables
ENV INPUT_DIR="/input"
ENV OUTPUT_DIR="/output"
# need to be mounted when running the container
ENV APP_DIR="/app"
ENV DINAMICA_EGO_CLI="/app/squashfs-root/usr/bin/DinamicaConsole"

# Create folders
RUN mkdir -p $INPUT_DIR $OUTPUT_DIR

# Install Micromamba
RUN apt update && apt install -y curl bzip2
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN eval "$(./bin/micromamba shell hook -s posix)"
RUN micromamba shell init -s bash -p ~/micromamba

# Create and activate environment - land_use
COPY 10_land_use_docker.sh /app/10_land_use_docker.sh
RUN chmod +x /app/10_land_use_docker.sh
COPY 10_land_use_requirements.txt /app/10_land_use_requirements.txt
RUN eval /app/10_land_use_docker.sh

# Download and Unpack Dinamica EGO 7 AppImage
# https://dinamicaego.com/nui_download/1792/
WORKDIR /app
RUN wget https://dinamicaego.com/nui_download/1792/ -O DinamicaEGO-740-Ubuntu-LTS.AppImage
RUN chmod +x DinamicaEGO-740-Ubuntu-LTS.AppImage
RUN ./DinamicaEGO-740-Ubuntu-LTS.AppImage --appimage-extract

# Define entrypoint for the DinamicaConsole in the AppImage, within Micromamba environment
ENTRYPOINT ["bin/bash", "-c", "source ~/micromamba/etc/profile.d/mamba.sh && source ~/micromamba/etc/profile.d/conda.sh && conda activate land_use && $DINAMICA_EGO_CLI"]

# Call with
# docker run -v /path/to/input:/input -v /path/to/output:/output dinamica_console:0.1 /path/to/model/model_file

# Delete docker container
# docker rm $(docker ps -a -q)