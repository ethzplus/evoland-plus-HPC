#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH --mem-per-cpu=10240
#SBATCH --tmp=8000            # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/10_land_use-%j.out"
#SBATCH --error="logs/10_land_use-%j.err"

echo "Current working directory: $(pwd)"
source /cluster/project/eawag/p01002/Future-EI/src/bash_common.sh
source /cluster/project/eawag/p01002/Future-EI/src/steps/10_LULCC/10_land_use.sh

# With tail -f logs/10_land_use-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/10_land_use-*.out logs/10_land_use-*.err