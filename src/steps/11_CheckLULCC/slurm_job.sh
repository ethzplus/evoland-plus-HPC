#!/bin/bash
#SBATCH --job-name="11_check_lulcc"
#SBATCH -n 1                  # Number of cores requested
#SBATCH --cpus-per-task=25    # Number of CPUs per task
#SBATCH --time=4:00:00        # Runtime
#SBATCH --mem-per-cpu=4G      # Memory per cpu in GB (see also --mem)
#SBATCH --tmp=2G              # https://scicomp.ethz.ch/wiki/Using_local_scratch
#SBATCH --output="logs/11_check_lulcc-%j.out"
#SBATCH --error="logs/11_check_lulcc-%j.err"
#SBATCH --mail-type=NONE       # Mail events (NONE, BEGIN, END, FAIL, ALL)

echo "Current working directory: $(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../de_activate.sh" check_lulc 1
echo "FUTURE_EI_CONFIG_FILE: $FUTURE_EI_CONFIG_FILE"
Rscript "$SCRIPT_DIR/11_CheckLULCC.R"

# With tail -f logs/11_check_lulcc-*.out you can follow the progress of the job
# To follow both .out and .err files, use tail -f logs/11_check_lulcc-*.out logs/11_check_lulcc-*.err

# Resources example: For 5 configurations, Switzerland map, sequential
# >> seff 3154855
# Nodes: 1
# Cores per node: 8
# CPU Utilized: 00:12:57
# CPU Efficiency: 11.96% of 01:48:16 core-walltime
# Job Wall-clock time: 00:13:32
# Memory Utilized: 3.38 GB
# Memory Efficiency: 42.21% of 8.00 GB

# Resources example: For 5 configurations, Switzerland map, concurrent
# >> seff 3158966
# Nodes: 1
# Cores per node: 2
# CPU Utilized: 00:08:51
# CPU Efficiency: 79.73% of 00:11:06 core-walltime
# Job Wall-clock time: 00:05:33
# Memory Utilized: 5.15 GB
# Memory Efficiency: 64.41% of 8.00 GB
# => S = 12:57/(2*8:51) = 0.75 efficient (should be higher with some more cores)