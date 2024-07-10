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
mark_finished_lulcc
log info "Got control table subset at $LULCC_M_SIM_CONTROL_TABLE"
log debug "Table: $(cat "$LULCC_M_SIM_CONTROL_TABLE")"

# Start Dinamica at different times, to avoid startup conflicts
# sleep for (SLURM_ARRAY_TASK_ID mod 50) seconds
wait_time=$(( SLURM_ARRAY_TASK_ID % 50 ))
log info "Sleeping for $wait_time seconds"
sleep $wait_time

log info "Starting LULCC step"
source "$FUTURE_EI_BASE_DIR"/src/steps/10_LULCC/10_land_use.sh

log info "Starting NCPs step"
source "$FUTURE_EI_BASE_DIR"/src/steps/40_NCPs/40_NCPs.sh

log info "Merge control table subset back to full control table"
merge_control_table


# With tail -f logs/10_40_combined_array_-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/10_40_combined_array_-*.out logs/10_40_combined_array_-*.err

# Resources example: For 1 configuration (part of array), Switzerland map
# >>> seff 62996867
# Nodes: 1
# Cores per node: 2
# CPU Utilized: 06:23:06
# CPU Efficiency: 75.29% of 08:28:50 core-walltime
# Job Wall-clock time: 04:14:25
# Memory Utilized: 7.56 GB
# Memory Efficiency: 75.58% of 10.00 GB

# For 1080 configurations and 4:30h max runtime this would result in 203 days
# of concurrent runtime. (Better assign a longer runtime to the job to avoid
# job termination due to time limit.)
#
# When using the array job, we can decide how many configurations we want to run
# in parallel. Slurm has a limit, which can be set by the cluster. To find out
# the limit, use `scontrol show config | grep MaxArraySize` and consider the
# resources available to the user.
# In our example, using 40 configurations in parallel, the runtime would be
# reduced to just above 5 days.