#!/bin/bash
#SBATCH --job-name="50_Summarisation"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=16    # Number of CPUs per task
#SBATCH --time=7-00:00:00     # Runtime in D-HH:MM:SS
#SBATCH --mem-per-cpu=4G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/50_summarisation-%j.out"
#SBATCH --error="logs/50_summarisation-%j.err"
#SBATCH --mail-type=NONE       # Mail events (NONE, BEGIN, END, FAIL, ALL)

echo "Current working directory: $(pwd)"
source "/cluster/project/eawag/p01002/Future-EI/src/de_activate.sh" ncp_summarisation 1
echo "FUTURE_EI_CONFIG_FILE: $FUTURE_EI_CONFIG_FILE"
Rscript /cluster/project/eawag/p01002/Future-EI/src/steps/50_Summarisation/50_Summarisation.R

# With tail -f logs/50_Summarisation-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/50_summarisation-*.out logs/50_summarisation-*.err

# Resources example: 