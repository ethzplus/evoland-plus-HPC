#!/bin/bash
#SBATCH --job-name="10_40_combined_array"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=4     # Number of CPUs per task
#SBATCH --time=0:20:00        # Runtime
#SBATCH --mem-per-cpu=4G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/10_40_combined_array-%j.out"
#SBATCH --error="logs/10_40_combined_array-%j.err"
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
## Array job
#SBATCH --array=1-2           # Number of array jobs - step size needs to be 1

echo "Current working directory: $(pwd)"
# $FUTURE_EI_BASE_DIR must be passed defined for array job, bash_common.sh
if [ -z "$FUTURE_EI_BASE_DIR" ]; then
  echo "FUTURE_EI_BASE_DIR is not defined. Exiting."
  exit 1
fi
source "$FUTURE_EI_BASE_DIR"/src/bash_common.sh

# Tell Dinamica max number of cores available
export DINAMICA_EGO_7_MAX_DETECTED_CORES
DINAMICA_EGO_7_MAX_DETECTED_CORES=$SLURM_CPUS_PER_TASK

mkdir -p "$FUTURE_EI_OUTPUT_DIR"

log info "Starting array job for LULCC and NCPs"
log info "Array job ID: $SLURM_ARRAY_JOB_ID"
log info "Array task ID: $SLURM_ARRAY_TASK_ID"
log info "Array count: $SLURM_ARRAY_TASK_COUNT"
# assure that the step size is 1
if [ "$SLURM_ARRAY_TASK_STEP" -ne 1 ]; then
  log error "Step size for array job needs to be 1. Got $SLURM_ARRAY_TASK_STEP."
  exit 1
fi

log debug "Current LULCC_M_SIM_CONTROL_TABLE: $LULCC_M_SIM_CONTROL_TABLE"
mkdir -p "$TMPDIR"

# Split the control table
source "$FUTURE_EI_BASE_DIR"/src/control_table_control.sh
split_control_table
log info "Got control table subset at $LULCC_M_SIM_CONTROL_TABLE"
log debug "Table: $(cat $LULCC_M_SIM_CONTROL_TABLE)"

#log info "Starting LULCC step"
#source "$FUTURE_EI_BASE_DIR"/src/steps/10_LULCC/10_land_use.sh

log info "Starting NCPs step"
source "$FUTURE_EI_BASE_DIR"/src/steps/40_NCPs/40_NCPs.sh

log info "Merge control table subset back to full control table"
merge_control_table


# With tail -f logs/10_40_combined_array_-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/10_40_combined_array_-*.out logs/10_40_combined_array_-*.err
