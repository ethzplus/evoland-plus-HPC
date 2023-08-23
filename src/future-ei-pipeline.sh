# Future EI Pipeline
# ------------------

# Variables
# ---------
tmp_dir=/tmp/$USER/$SLURM_JOB_ID  # Temporary directory
keep_tmp=0                        # Keep temporary files (0: no, 1: yes) for debugging

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
# Max time: days-hours:minutes:seconds
# SBATCH --time=0-00:02:00
## Array job
## SBATCH --array=1-4%2

# Log file
# SBATCH --output=logs/future-ei-%j.out
# Error file
# SBATCH --error=logs/future-ei-%j.err
# Email notifications
# SBATCH --mail-type=ALL,TIME_LIMIT_80
# SBATCH --mail-user=<email>

# Change directory before running
# SBATCH --chdir=~/Future-EI/src

# Find Further Environment Variables here
# https://slurm.schedmd.com/sbatch.html#SECTION_INPUT-ENVIRONMENT-VARIABLES

# ______________________________________________________________________________________

srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE hostname

# Preparation
# -----------
# Make temporary directories - one for every node
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE mkdir -p $tmp_dir
# Preparation script
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/00_Preparation.mpi

# ______________________________________________________________________________________


# LULCC
# -----
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/01_LULCC.mpi

# Focal LULC Preparation
# ----------------------
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/20_FocalPrep.mpi

# N-SDM Biodiversity
# ------------------
# all $SLURM_JOB_NUM_NODES x $SLURM_TASKS_PER_NODE tasks will be used
srun steps/30_N-SDM.mpi

# NCPs
# ----
srun steps/40_NCPs.mpi

# Species Aggregation
# -------------------
srun steps/31_SpeciesAgg.mpi &

# NCP Clustering
# --------------
srun steps/41_NCP_Clustering.mpi &

# last two steps can be run in parallel (caused by &)

# ______________________________________________________________________________________

# Clean up
# --------
# Save results to output directory
srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE steps/99_Cleanup.mpi
# Remove temporary directories
if [ $keep_tmp -eq 0 ]; then
    srun --nodes=$SLURM_JOB_NUM_NODES --tasks-per-node=$SLURM_TASKS_PER_NODE rm -rf $tmp_dir
fi

wait  # Wait for all srun commands to finish before marking as done