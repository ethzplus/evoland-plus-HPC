#!/bin/bash
#SBATCH --job-name="50_analysis"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=12    # Number of CPUs per task
#SBATCH --time=7-00:00:00     # Runtime in D-HH:MM:SS
#SBATCH --mem-per-cpu=6G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/50_analysis-%j.out"
#SBATCH --error="logs/50_analysis-%j.err"
#SBATCH --mail-type=NONE       # Mail events (NONE, BEGIN, END, FAIL, ALL)

echo "Current working directory: $(pwd)"
source "/cluster/project/eawag/p01002/Future-EI/src/de_activate.sh" ncp_summarisation 1
echo "FUTURE_EI_CONFIG_FILE: $FUTURE_EI_CONFIG_FILE"
Rscript /cluster/project/eawag/p01002/Future-EI/src/steps/50_Analysis/51_Summarisation.R
#Rscript /cluster/project/eawag/p01002/Future-EI/src/steps/50_Analysis/52_Clustering.R

# With tail -f logs/50_analysis-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/50_analysis-*.out logs/50_analysis-*.err

# Resources example: 