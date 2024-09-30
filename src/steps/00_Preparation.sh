#!/bin/bash
# Preparation Job

# Use scratch for output folders
log debug "Using scratch for output folders:"
output_dirs=(
    LULCC_CH_OUTPUT_BASE_DIR
    CHECK_LULCC_OUTPUT_DIR
    NCP_OUTPUT_BASE_DIR
    NSDM_OUTPUT_BASE_DIR
    FOCAL_OUTPUT_BASE_DIR
    SPAGG_OUTPUT_BASE_DIR
)
for dir in "${output_dirs[@]}"; do
    declare "$dir"="$FUTURE_EI_OUTPUT_DIR${!dir}"
    log debug "$dir: ${!dir}"
done
scratch_dirs=(
    NCP_SCRATCH_BASE_DIR
)
for dir in "${scratch_dirs[@]}"; do
    declare "$dir"="$TMPDIR${!dir}"
    log debug "$dir: ${!dir}"
done
