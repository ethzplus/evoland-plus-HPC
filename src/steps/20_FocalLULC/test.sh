#!/bin/bash
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=20     # Number of CPUs per task
#SBATCH --time=1:00:00        # Runtime
#SBATCH --mem-per-cpu=4G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/20_focal_statistics-%j.out"
#SBATCH --error="logs/20_focal_statistics-%j.err"
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)

echo "Current working directory: $(pwd)"
source "/cluster/project/eawag/p01002/Future-EI/src/de_activate.sh" focal_lulc 1
echo "FUTURE_EI_CONFIG_FILE: $FUTURE_EI_CONFIG_FILE"
Rscript /cluster/project/eawag/p01002/Future-EI/src/steps/20_FocalLULC/20_focal_statistics.R

# With tail -f logs/20_focal_statistics-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/20_focal_statistics-*.out logs/20_focal_statistics-*.err