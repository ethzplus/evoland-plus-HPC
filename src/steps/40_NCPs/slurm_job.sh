#!/bin/bash
#SBATCH --job-name="40_NCPs_multiple"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=6     # Number of CPUs per task
#SBATCH --time=1:00:00        # Runtime
#SBATCH --mem-per-cpu=3G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=12G             # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/40_NCPs_multiple-%j.out"
#SBATCH --error="logs/40_NCPs_multiple-%j.err"
#SBATCH --mail-type=NONE      # Mail events (NONE, BEGIN, END, FAIL, ALL)

# NCP Slurm Job
# -------------
# This script runs all NCP models for all scenarios and years in the
# scenario table.
# The parallelized LULCC and NCP script can be found one folder up
# '10_40_combined_array_job.sh'.

echo "Current working directory: $(pwd)"
source /cluster/project/eawag/p01002/Future-EI/src/bash_common.sh

source /cluster/project/eawag/p01002/Future-EI/src/steps/40_NCPs/40_NCPs.sh

# With tail -f logs/40_NCPs_multiple-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/40_NCPs_multiple-*.out logs/40_NCPs_multiple-*.err

# Resources example: For 5 configurations, Switzerland map
# >>> seff