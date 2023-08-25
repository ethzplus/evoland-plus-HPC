#!/bin/bash
# ------------------
# Future EI Pipeline
# ------------------

source logging.sh

# Script Variables
# ----------------
TMP_DIR=$TMPDIR/$SLURM_JOB_ID        # Temporary directory
KEEP_TMP=0                           # Keep temporary files (0: no, 1: yes) for debug
LOG_DIR=~/logs                       # Log directory

# Slurm
# -----
# Request 3 nodes with 2 task each -> 6 CPUs
# SBATCH --nodes=3
# SBATCH --tasks-per-node=2
## Number of tasks
## SBATCH --ntasks=6
# Memory per node
# SBATCH --mem=1G
## Memory per CPU
## SBATCH --mem-per-cpu=500M
# Scratch space
# SBATCH --tmp=1G
# Max time: days-hours:minutes:seconds
# SBATCH --time=0-00:02:00
## Array job
## SBATCH --array=1-4%2

# Log file
# SBATCH --output=$LOG_DIR/future-ei-%j.out
# Error file
# SBATCH --error=$LOG_DIR/future-ei-%j.err
# Email notifications
# SBATCH --mail-type=ALL,TIME_LIMIT_90
# SBATCH --mail-user=<email>

# Change directory before running
# SBATCH --chdir=~/Future-EI/src

# Find Further Environment Variables here
# https://slurm.schedmd.com/sbatch.html#SECTION_INPUT-ENVIRONMENT-VARIABLES

# ______________________________________________________________________________________

log info "Starting Future EI Pipeline"
log info "Job ID: $SLURM_JOB_ID"
log info "Job Name: $SLURM_JOB_NAME"
if [ ! -z "$SLURM_ARRAY_JOB_ID" ]; then
    log info "Job Array ID: $SLURM_ARRAY_JOB_ID"
    log info "Job Array Index: $SLURM_ARRAY_TASK_ID"
fi
log info "Number of Nodes: $SLURM_JOB_NUM_NODES"
log info "Number of Tasks per Node: $SLURM_TASKS_PER_NODE"
log info "Temporary Directory: $TMP_DIR"

srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE hostname

# Preparation
# -----------
# Make temporary directories - one for every node
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE mkdir -p $TMP_DIR
# Preparation script
log info "Running Preparation - 00_Preparation.sh"
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/00_Preparation.sh

# ______________________________________________________________________________________


# 01_LULCC - LULCC
# ----------------
log info "Running LULCC - 01_LULCC.sh"
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/01_LULCC.sh

# 20_FocalPrep - Focal LULC Preparation
# -------------------------------------
log info "Running Focal LULC Preparation - 20_FocalPrep.sh"
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/20_FocalPrep.sh

# 30_N-SDM - N-SDM Biodiversity
# -----------------------------
# all $SLURM_JOB_NUM_NODES x $SLURM_TASKS_PER_NODE tasks will be used
log info "Running N-SDM Biodiversity - 30_N-SDM.sh - may take a while"
srun steps/30_N-SDM.sh

# 40_NCPs - NCPs
# --------------
log info "Running NCPs - 40_NCPs.sh"
srun steps/40_NCPs.sh

# 31_SpeciesAgg - Species Aggregation
# -----------------------------------
log info "Running Species Aggregation - 31_SpeciesAgg.sh"
srun steps/31_SpeciesAgg.sh &

# 41_NCP_Clustering - NCP Clustering
# ----------------------------------
log info "Running NCP Clustering - 41_NCP_Clustering.sh"
srun steps/41_NCP_Clustering.sh &

# last two steps can be run in parallel (caused by &)

# ______________________________________________________________________________________

# 99_Cleanup - Clean up
# ---------------------
# Save results to output directory
log info "Running Cleanup - 99_Cleanup.sh"
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/99_Cleanup.sh
# Remove temporary directories
if [ $KEEP_TMP -eq 0 ]; then
    srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE rm -rf $TMP_DIR
else
    log warning "Keeping temporary files in $TMP_DIR. Non-permanent if in local scratch."
fi

wait  # Wait for all srun commands to finish before marking as done