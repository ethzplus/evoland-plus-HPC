# Docker image the LULCC use

# Use Ubuntu 22.10
FROM ubuntu:23.04
LABEL authors="Carlson Büth" \
      version="0.1" \
      description="Docker image for the LULCC (land use land cover change) model"

# Environment variables
ENV INPUT_DIR="/input"
ENV OUTPUT_DIR="/output"
ENV MODEL_DIR="/model"
# need to be mounted when running the container
ENV APP_DIR="/app"
ENV DINAMICA_EGO_CLI="$APP_DIR/squashfs-root/usr/bin/DinamicaConsole"
# Add shared libraries to the library path
ENV LD_LIBRARY_PATH="$APP_DIR/squashfs-root/usr/lib/:$LD_LIBRARY_PATH"
ENV DINAMICA_EGO_7_INSTALLATION_DIRECTORY="$APP_DIR/squashfs-root/usr/bin"
ENV REGISTRY_FILE="/root/.dinamica_ego_7.conf"

# Create folders
RUN mkdir -p $INPUT_DIR $OUTPUT_DIR $MODEL_DIR

# Install Micromamba
RUN apt update && apt install -y curl bzip2
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN eval "$(./bin/micromamba shell hook -s posix)"
RUN micromamba shell init -s bash -p ~/micromamba

# Create and activate environment - land_use
COPY 10_land_use_docker.sh $APP_DIR/10_land_use_docker.sh
RUN chmod +x $APP_DIR/10_land_use_docker.sh
COPY 10_land_use_requirements.txt $APP_DIR/10_land_use_requirements.txt
RUN eval $APP_DIR/10_land_use_docker.sh

# Download and Unpack Dinamica EGO 7 AppImage
# https://dinamicaego.com/nui_download/1792/
WORKDIR $APP_DIR
RUN curl -Ls https://dinamicaego.com/nui_download/1792/ -o DinamicaEGO-740-Ubuntu-LTS.AppImage
RUN chmod +x DinamicaEGO-740-Ubuntu-LTS.AppImage
RUN ./DinamicaEGO-740-Ubuntu-LTS.AppImage --appimage-extract

WORKDIR $MODEL_DIR
COPY .dinamica_ego_7.conf $REGISTRY_FILE

# Define entrypoint for the DinamicaConsole in the AppImage, within Micromamba environment
#ENTRYPOINT ["/bin/bash", "-c", "eval \"$(micromamba shell hook --shell bash)\" && micromamba activate land_use && exec $DINAMICA_EGO_CLI $APP_DIR/LULCC_Dinamica_models/LULCC_CH.ego \"$@\"", "--"]
#ENTRYPOINT ["/bin/bash", "-c", "eval \"$(micromamba shell hook --shell bash)\" && micromamba activate land_use && exec $DINAMICA_EGO_CLI \"$@\" $MODEL_DIR/$MODEL_SCRIPT", "--"]
#Rscript Scripts/LULCC_CH_master.R
ENTRYPOINT ["/bin/bash", "-c", "eval \"$(micromamba shell hook --shell bash)\" && micromamba activate land_use && exec Rscript $MODEL_DIR/Scripts/LULCC_CH_master.R", "--"]


# Build
# docker build -t lulcc:0.1 .

# Call with
# docker run -v /path/to/model:/model -e MODEL_SCRIPT="/app/model_file"  -v /path/to/input:/input -v /path/to/output:/output lulcc:0.1

# Delete docker container
# docker rm $(docker ps -a -q)