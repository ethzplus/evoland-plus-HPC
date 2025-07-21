#!/bin/bash
# NCP step for the batch job src/10_40_combined_array_job.sh
# For each simulation and year, run the NCP models
# Uses run_all_ncps.sh <NCP_RUN_SCENARIO_ID> <NCP_RUN_YEAR> <NCP_RUN_INPUT_DIR> <NCP_RUN_OUTPUT_DIR> <NCP_RUN_SCRATCH_DIR>
# If called separately for testing, source the bash_common.sh script first

# $FUTURE_EI_BASE_DIR must be passed, as this script is used in array jobs
if [ -z "$FUTURE_EI_BASE_DIR" ]; then
  echo "FUTURE_EI_BASE_DIR is not defined. Exiting."
  exit 1
fi
source "$FUTURE_EI_BASE_DIR"/src/bash_common.sh
log info "Starting NCPs step"
# Scenario IDs are taken from the simulation control table,
# which is found at $LULCC_M_SIM_CONTROL_TABLE (might be subset for array job)
log debug "Current LULCC_CH_OUTPUT_BASE_DIR: $LULCC_CH_OUTPUT_BASE_DIR"
log debug "Current LULCC_M_SIM_CONTROL_TABLE: $LULCC_M_SIM_CONTROL_TABLE"
log debug "Current NCP_OUTPUT_BASE_DIR: $NCP_OUTPUT_BASE_DIR"
# Get the simulation IDs from the control table
# Simulation_num.,Scenario_ID.string,Simulation_ID.string,Model_mode.string,
# Scenario_start.real,Scenario_end.real,Step_length.real,...
simulation_ids=$(tail -n +2 "$LULCC_M_SIM_CONTROL_TABLE" | cut -d, -f3)
# tail  -n +2: start from the second line, skipping the header
# cut -d, -f3: get the third field, which is the simulation ID
log info "Found $(echo "$simulation_ids" | wc -w) simulation IDs in the current control table"

# Create the output directories if they do not exist
mkdir -p "$NCP_OUTPUT_BASE_DIR"
mkdir -p "$NCP_SCRATCH_BASE_DIR"

# Loop over simulation IDs
for simulation_id in $simulation_ids; do
  log info "Running NCPs for Simulation ID $simulation_id"
  # Instead of using the year from the control table, the available lulc output
  # layers are used $LULCC_CH_OUTPUT_BASE_DIR:
  # {simulation_id}/simulated_LULC_simID_{simulation_id}_year_{year}.tif
  # shellcheck disable=SC2010
  years=$(ls "$LULCC_CH_OUTPUT_BASE_DIR/$simulation_id" | \
          grep -oP "simulated_LULC_simID_${simulation_id}_year_\K\d+(?=\.tif$)"\
          | sort -n)
  # ls      : list the files in the simulation output directory
  # grep -oP: extract the year from the file name
  # sort  -n: sort the years numerically
  log info "With $(echo "$years" | wc -w) years: [$(echo "$years" | tr '\n' ' ')]"
  # Loop over years
  for year in $years; do
    log info "Running NCPs for simulation ID $simulation_id and year $year"
    # Run the NCP models
    source "$FUTURE_EI_BASE_DIR"/src/steps/40_NCPs/NCP_models/run_all_ncps.sh \
           "$simulation_id" "$year" \
           "$LULCC_CH_OUTPUT_BASE_DIR" \
           "$NCP_OUTPUT_BASE_DIR" \
           "$NCP_SCRATCH_BASE_DIR"
  done
done
log info "Finished NCPs step"
