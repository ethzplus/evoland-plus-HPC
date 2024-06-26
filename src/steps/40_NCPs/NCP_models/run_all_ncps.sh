#!/bin/bash
# Execute all NCPs with one LULC map/year
# Usage: bash run_all_ncps.sh <NCP_RUN_SCENARIO_ID> <NCP_RUN_YEAR> <NCP_RUN_INPUT_DIR> <NCP_RUN_OUTPUT_DIR> <NCP_RUN_SCRATCH_DIR>
# Example: bash run_all_ncps.sh 1 2015 /path/to/input_dir /path/to/output_dir /path/to/scratch_dir
# ------------------------------------------------------------------------------

# Script Variables
# ----------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../../bash_common.sh"
# Activate the conda environment with R and python packages
log debug "Activating conda environment: ncps"
source "$SCRIPT_DIR/../../../de_activate.sh" ncps 1
# Check if the conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "ncps" ]; then
    log error "The conda environment is not activated: ncps"
    exit 1
fi

source "$SCRIPT_DIR"/script_utils.sh


# ------------------------------------------------------------------------------

# Check input arguments
if [ "$#" -ne 5 ]; then
    log error "Illegal number of parameters ($#)"
    log error "Usage: bash run_all_ncps.sh <NCP_RUN_SCENARIO_ID> <year> <NCP_RUN_INPUT_DIR> <NCP_RUN_OUTPUT_DIR> <NCP_RUN_SCRATCH_DIR>"
    log error "Received: $*"
    exit 1
fi

# Assign input arguments as environmental variables
export NCP_RUN_SCENARIO_ID
export NCP_RUN_YEAR
export NCP_RUN_INPUT_DIR
export NCP_RUN_OUTPUT_DIR
export NCP_RUN_SCRATCH_DIR
NCP_RUN_SCENARIO_ID=$1
NCP_RUN_YEAR=$2
NCP_RUN_INPUT_DIR=$3
NCP_RUN_OUTPUT_DIR=$4
NCP_RUN_SCRATCH_DIR=$5

# Check if the input scenario ID is a number
if ! [[ "$NCP_RUN_SCENARIO_ID" =~ ^[0-9]+$ ]]; then
    log error "The scenario ID must be a number: $NCP_RUN_SCENARIO_ID"
    exit 1
fi

# Check if the year is a number
if ! [[ "$NCP_RUN_YEAR" =~ ^[0-9]+$ ]]; then
    log error "The year must be a number: $NCP_RUN_YEAR"
    exit 1
fi

# Check if the input directory exists
if [ ! -d "$NCP_RUN_INPUT_DIR" ]; then
    log error "The input directory does not exist: $NCP_RUN_INPUT_DIR"
    exit 1
fi

# Check if the output directory exists
if [ ! -d "$NCP_RUN_OUTPUT_DIR" ]; then
    log error "The output directory does not exist: $NCP_RUN_OUTPUT_DIR"
    exit 1
fi

# Check if the scratch directory exists
if [ ! -d "$NCP_RUN_SCRATCH_DIR" ]; then
    log error "The scratch directory does not exist: $NCP_RUN_SCRATCH_DIR"
    exit 1
fi

# ------------------------------------------------------------------------------

# Set and copy NCP_PARAMS_YML to NCP_RUN_OUTPUT_DIR
new_ncp_params_path="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/run_ncp_params_$NCP_RUN_YEAR.yml"
python "$SCRIPT_DIR/load_params.py" "$new_ncp_params_path"
export NCP_PARAMS_YML
NCP_PARAMS_YML_BASE=$NCP_PARAMS_YML
NCP_PARAMS_YML=$new_ncp_params_path


# Models
log info "Running all NCPs: Scenario ID=$NCP_RUN_SCENARIO_ID, Year=$NCP_RUN_YEAR, Input dir=$NCP_RUN_INPUT_DIR, Output dir=$NCP_RUN_OUTPUT_DIR, Scratch dir=$NCP_RUN_SCRATCH_DIR"


## CAR - Regulation of climate
## Indicator: Carbon sequestration
log info "Running CAR - Carbon sequestration"
run_scripts "$SCRIPT_DIR/CAR/1_CAR_S_CH.R" \
            "$SCRIPT_DIR/CAR/2_CAR_S_CH.py" \
            "$SCRIPT_DIR/CAR/3_CAR_S_CH.R"

## FF - Food and feed
## Indicator: Crop production potential (`ecocrop`)
log info "Running FF - Crop production potential"
# needs FF/data_preparation/0_monthly_average_dataprep.R from prepare_ncps.sh
run_scripts "$SCRIPT_DIR/FF/1_FF_aggregate_masking.R"

## HAB - Habitat creation and maintenance
## Indicator: Habitat quality index
log info "Running HAB - Habitat quality index"
run_scripts "$SCRIPT_DIR/HAB/0_HAB_threat_layers_generation.R" \
            "$SCRIPT_DIR/HAB/1_HAB_S_CH.py"

## NDR - Nutrient Delivery Ratio
## Indicator: Annual nutrient retention by vegetation
log info "Running NDR - Annual nutrient retention by vegetation"
run_scripts "$SCRIPT_DIR/NDR/1_NDR_S_CH.py"

## POL - Pollination and dispersal of seeds
## Indicator: Habitat abundance for pollinators
log info "Running POL - Habitat abundance for pollinators"
run_scripts "$SCRIPT_DIR/POL/1_POL_S_CH.py" \
            "$SCRIPT_DIR/POL/2_POL_S_CH_aggregating.R"

## SDR - Formation, protection and decontamination of soils
## Indicator: Erosion control by sediment retention
log info "Running SDR - Erosion control by sediment retention"
run_scripts "$SCRIPT_DIR/SDR/1_SDR_S_CH.py"

## WY - Regulation of freshwater quantity, location and timing
## Indicator: Annual water yield
log info "Running WY - Annual water yield"
run_scripts "$SCRIPT_DIR/WY/1_WY_S_CH.py"

# ------------------------------------------------------------------------------

log info "All NCPs have been run"

NCP_PARAMS_YML=$NCP_PARAMS_YML_BASE