#!/bin/bash

#SBATCH --nodes=1            # Number of nodes requested
#SBATCH --tasks-per-node=4   # Number of cores per node
#SBATCH --cpus-per-task=4    # Number of CPUs per task
# Memory per node
#SBATCH --mem-per-cpu=4G     # Memory per processor
#SBATCH --time=6:00:00        # Runtime
# Scratch space
#SBATCH --tmp=12G            # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/10_land_use-%j.out"
#SBATCH --error="logs/10_land_use-%j.err"
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)

# [LULCC_START_ROW, LULCC_END_ROW[ is the range of rows in the
# Simulation_control.csv file that will be processed by this job
export LULCC_START_ROW
export LULCC_END_ROW
LULCC_START_ROW=3
LULCC_END_ROW=5

# Tell Dinamica max number of cores available
export DINAMICA_EGO_7_MAX_DETECTED_CORES
DINAMICA_EGO_7_MAX_DETECTED_CORES=$SLURM_CPUS_PER_TASK

echo "Current working directory: $(pwd)"
source /cluster/project/eawag/p01002/Future-EI/src/bash_common.sh
source /cluster/project/eawag/p01002/Future-EI/src/steps/10_LULCC/10_land_use.sh

mkdir -p "$FUTURE_EI_OUTPUT_DIR"
cp -r "$LULCC_CH_OUTPUT_BASE_DIR" "$FUTURE_EI_OUTPUT_DIR/land_use_test"

# With tail -f logs/10_land_use-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/10_land_use-*.out logs/10_land_use-*.err

# Resources example: For 5 configurations, Switzerland map
# >>> seff 51490116
# Nodes: 1
# Cores per node: 4
# CPU Utilized: 10:53:37
# CPU Efficiency: 49.48% of 22:01:00 core-walltime
# Job Wall-clock time: 05:30:15 -> per configuration: 1:06:03
# Memory Utilized: 8.16 GB
# Memory Efficiency: 50.99% of 16.00 GB