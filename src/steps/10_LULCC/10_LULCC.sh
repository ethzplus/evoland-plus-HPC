#!/bin/bash
# Land Use Land Cover Change Job

# Needs docker image lulcc - build or pull it first
# docker build -t lulcc:0.1 .
# or
# docker pull lulcc/lulcc:0.1
docker run -v "$LULCC_CH_HPC_DIR":/model lulcc:0.1