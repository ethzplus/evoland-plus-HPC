#!/bin/bash
# ------------------
# Future EI Pipeline
# ------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "$SCRIPT_DIR"/bash_common.sh

# Slurm
# -----
# Request 3 nodes with 2 task each -> 6 CPUs
# SBATCH --nodes=3
# SBATCH --tasks-per-node=2
# Memory per node
# SBATCH --mem=1G
# Scratch space
# SBATCH --tmp=1G
# Max time: days-hours:minutes:seconds
# SBATCH --time=0-00:02:00

# Log file
# SBATCH --output=~/logs/future-ei-%j.out
# Error file
# SBATCH --error=~/logs/future-ei-%j.err
# Email notifications
# SBATCH --mail-type=ALL,TIME_LIMIT_90
# SBATCH --mail-user=<email>

## Change directory before running
## SBATCH --chdir=~/Future-EI/src

# Find Further Environment Variables here
# https://slurm.schedmd.com/sbatch.html#SECTION_INPUT-ENVIRONMENT-VARIABLES

# ______________________________________________________________________________________

log info "Starting Future EI Pipeline"
log info "Job ID: $SLURM_JOB_ID"
log info "Job Name: $SLURM_JOB_NAME"
log info "Number of Nodes: $SLURM_JOB_NUM_NODES"
log info "Number of Tasks per Node: $SLURM_TASKS_PER_NODE"
log info "Temporary Directory: $TMPDIR"

# ______________________________________________________________________________________


srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" hostname

# Preparation
# -----------
# Make temporary directories - one for every node
srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" mkdir -p "$TMPDIR"
## Preparation script
#log info "Running Preparation - 00_Preparation.sh"
#srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" steps/00_Preparation.sh


# only if array job
if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
    source "$SCRIPT_DIR"/control_table_control.sh
    # Split control table - get rows corresponding to this job array
    split_control_table
    log info "Got control table subset at $LULCC_M_SIM_CONTROL_TABLE"
else
    # set LULCC_START_ROW and LULCC_END_ROW to 2 and max row number
    LULCC_START_ROW=2
    LULCC_END_ROW=$(($(grep -c -v '^$' "$LULCC_M_SIM_CONTROL_TABLE")+1))
    export LULCC_START_ROW
    export LULCC_END_ROW
    log info "Running all rows of $LULCC_M_SIM_CONTROL_TABLE"
fi


# ______________________________________________________________________________________


# 10_LULCC - LULCC
# ----------------
mkdir -p "$LULCC_CH_OUTPUT_BASE_DIR"
log info "Running LULCC - 10_LULCC.sh"
srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" steps/10_LULCC/10_LULCC.sh

# Remove prediction probability maps if
# $LULCC_M_REMOVE_PRED_PROB_MAPS is 1 or True/TRUE
if [ "$LULCC_M_REMOVE_PRED_PROB_MAPS" = "1" ] || [ "$LULCC_M_REMOVE_PRED_PROB_MAPS" = "True" ] || [ "$LULCC_M_REMOVE_PRED_PROB_MAPS" = "TRUE" ]; then
    log info "Removing prediction probability maps"
    rm -rf "$LULCC_CH_HPC_DIR/Results/Pred_prob_maps"
fi


# only if array job
if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
    # Merge control table back together
    merge_control_table
fi

# Moving land use layers to N-SDM folder

# loop over rows - loop_id
for i in {LULCC_START_ROW..LULCC_END_ROW}; do
  log info "Running loop $i (parameters: $(sed -n "$i"p "$LULCC_M_SIM_CONTROL_TABLE"))"
  
  SIM_ID=$i
  LULCC_CH_OUTPUT_SIM_DIR=$LULCC_CH_OUTPUT_BASE_DIR$SIM_ID  # TODO: is SIM_ID needed?
  FOCAL_OUTPUT_SIM_DIR=$FOCAL_OUTPUT_BASE_DIR$SIM_ID
  log debug "Current LULCC_CH_OUTPUT_SIM_DIR: $LULCC_CH_OUTPUT_SIM_DIR"
  log debug "Current FOCAL_OUTPUT_SIM_DIR: $FOCAL_OUTPUT_SIM_DIR"

  # 20_FocalPrep - Focal LULC Preparation
  # -------------------------------------
  log info "Running Focal LULC Preparation - 20_FocalPrep.sh"
  srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" steps/20_FocalLULC/20_FocalPrep.sh
  # TODO: include transformation and reorganization routine + renaming and folder structure

  NSDM_OUTPUT_SIM_DIR=$NSDM_OUTPUT_BASE_DIR$SIM_ID
  log debug "Current NSDM_OUTPUT_SIM_DIR: $NSDM_OUTPUT_SIM_DIR"

  # 30_N-SDM - N-SDM Biodiversity
  # -----------------------------
  # all $SLURM_JOB_NUM_NODES x $SLURM_TASKS_PER_NODE tasks will be used
  log info "Running N-SDM Biodiversity - 30_N-SDM.sh - may take a while"
  srun steps/30_N-SDM.sh

  SPAGG_OUTPUT_SIM_DIR=$SPAGG_OUTPUT_BASE_DIR$SIM_ID
  log debug "Current SPAGG_OUTPUT_SIM_DIR: $SPAGG_OUTPUT_SIM_DIR"

  # 31_SpeciesAgg - Species Aggregation
  # -----------------------------------
  log info "Running Species Aggregation - 31_SpeciesAgg.sh"
  srun steps/31_SpeciesAgg.sh #&

  NCP_OUTPUT_SIM_DIR="$NCP_OUTPUT_BASE_DIR/$SIM_ID"
  log debug "Current NCP_OUTPUT_SIM_DIR: $NCP_OUTPUT_SIM_DIR"

  # 40_NCPs - NCPs
  # --------------
  log info "Running NCPs - 40_NCPs.sh"
  srun steps/40_NCPs.sh

  
  # 41_NCP_Clustering - NCP Clustering
  # ----------------------------------
  #log info "Running NCP Clustering - 41_NCP_Clustering.sh"
  #srun steps/41_NCP_Clustering.sh #& # for parallel execution

done

# ______________________________________________________________________________________

# 99_Cleanup - Clean up
# ---------------------
# Save results to output directory
log info "Running Cleanup - 99_Cleanup.sh"
#srun --nodes="$SLURM_JOB_NUM_NODES" --tasks-per-node="$SLURM_TASKS_PER_NODE" steps/99_Cleanup.sh

wait  # Wait for all srun commands to finish before marking as done