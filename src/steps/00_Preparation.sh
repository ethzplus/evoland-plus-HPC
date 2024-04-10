#!/bin/bash
# Preparation Job

# Use scratch for output folders
log debug "Using scratch for output folders:"
output_dirs=(
    LULCC_CH_OUTPUT_BASE_DIR
    NCP_OUTPUT_BASE_DIR
    NSDM_OUTPUT_BASE_DIR
    FOCAL_OUTPUT_BASE_DIR
    SPAGG_OUTPUT_BASE_DIR
)
for dir in "${output_dirs[@]}"; do
    declare "$dir"="$TMPDIR${!dir}"
    log debug "$dir: ${!dir}"
done