#!/bin/bash

#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=8     # Number of CPUs per task
#SBATCH --time=4:00:00        # Runtime
#SBATCH --mem-per-cpu=16G     # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/10_land_use-%j.out"
#SBATCH --error="logs/10_land_use-%j.err"
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)

echo "Current working directory: $(pwd)"
source /cluster/project/eawag/p01002/Future-EI/src/bash_common.sh
source /cluster/project/eawag/p01002/Future-EI/src/steps/10_LULCC/10_land_use.sh

# With tail -f logs/10_land_use-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/10_land_use-*.out logs/10_land_use-*.err