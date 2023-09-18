#!/bin/bash
# Land Use Land Cover Change Job

# Needs docker image lulcc
docker run -v $LULCC_CH_HPC_DIR:/model lulcc:0.1