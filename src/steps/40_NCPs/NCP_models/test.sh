#!/bin/bash
#SBATCH --job-name="40_NCPs_single"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=1     # Number of CPUs per task
#SBATCH --time=4:00:00        # Runtime
#SBATCH --mem-per-cpu=20G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/40_NCPs_single-%j.out"
#SBATCH --error="logs/40_NCPs_single-%j.err"
#SBATCH --mail-type=NONE      # Mail events (NONE, BEGIN, END, FAIL, ALL)

year=2020
scenario_id=1

echo "Current working directory: $(pwd)"
source /cluster/project/eawag/p01002/Future-EI/src/bash_common.sh
mkdir -p "$NCP_OUTPUT_BASE_DIR"
mkdir -p "$NCP_SCRATCH_BASE_DIR"

#source /cluster/project/eawag/p01002/Future-EI/src/steps/40_NCPs/NCP_models/prepare_ncps.sh
source /cluster/project/eawag/p01002/Future-EI/src/steps/40_NCPs/NCP_models/run_all_ncps.sh \
       $scenario_id $year \
       "$LULCC_CH_OUTPUT_BASE_DIR" \
       "$NCP_OUTPUT_BASE_DIR" \
       "$NCP_SCRATCH_BASE_DIR"


# With tail -f logs/40_NCPs_single-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/40_NCPs_single-*.out logs/20_focal_statistics-*.err
