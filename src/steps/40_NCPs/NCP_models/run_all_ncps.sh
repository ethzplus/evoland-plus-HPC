#!/bin/bash
# Execute all NCPs with automatic environment switching
# Usage: bash run_all_ncps.sh <NCP_RUN_SCENARIO_ID> <NCP_RUN_YEAR> <NCP_RUN_INPUT_DIR> <NCP_RUN_OUTPUT_DIR> <NCP_RUN_SCRATCH_DIR>

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Running script from: $SCRIPT_DIR"
source "$SCRIPT_DIR/../../../bash_common.sh"
source "$SCRIPT_DIR"/script_utils.sh

# Function to activate environment and check
activate_env() {
    local env_name=$1
    log debug "Activating conda environment: $env_name"
    source "$SCRIPT_DIR/../../../de_activate.sh" "$env_name" 1
    if [ "$CONDA_DEFAULT_ENV" != "$env_name" ]; then
        log error "Failed to activate conda environment: $env_name"
        exit 1
    fi
}

# Function to run scripts with appropriate environment
run_scripts_with_env() {
    local env_type=$1
    shift
    local scripts=("$@")
    
    case $env_type in
        "python")
            activate_env "ncps-python"
            ;;
        "r")
            activate_env "ncps-r"
            ;;
        *)
            log error "Unknown environment type: $env_type"
            exit 1
            ;;
    esac
    
    for script in "${scripts[@]}"; do
        log debug "Running script: $script"
        if [[ $script == *.py ]]; then
            python "$script"
        elif [[ $script == *.R ]]; then
            Rscript "$script"
        else
            log error "Unknown script type: $script"
            exit 1
        fi
    done
}

# Check input arguments (same as before)
if [ "$#" -ne 5 ]; then
    log error "Illegal number of parameters ($#)"
    log error "Usage: bash run_all_ncps.sh <NCP_RUN_SCENARIO_ID> <year> <NCP_RUN_INPUT_DIR> <NCP_RUN_OUTPUT_DIR> <NCP_RUN_SCRATCH_DIR>"
    log error "Received: $*"
    exit 1
fi

# Assign input arguments (same as before)
export NCP_RUN_SCENARIO_ID=$1
export NCP_RUN_YEAR=$2
export NCP_RUN_INPUT_DIR=$3
export NCP_RUN_OUTPUT_DIR=$4
export NCP_RUN_SCRATCH_DIR=$5

# Validation (same as before)
if ! [[ "$NCP_RUN_SCENARIO_ID" =~ ^[0-9]+$ ]]; then
    log error "The scenario ID must be a number: $NCP_RUN_SCENARIO_ID"
    exit 1
fi

if ! [[ "$NCP_RUN_YEAR" =~ ^[0-9]+$ ]]; then
    log error "The year must be a number: $NCP_RUN_YEAR"
    exit 1
fi

for dir in "$NCP_RUN_INPUT_DIR" "$NCP_RUN_OUTPUT_DIR" "$NCP_RUN_SCRATCH_DIR"; do
    if [ ! -d "$dir" ]; then
        log error "Directory does not exist: $dir"
        exit 1
    fi
done

# Set up parameters (using Python environment)
activate_env "ncps-python"
new_ncp_params_path="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/run_ncp_params_$NCP_RUN_YEAR.yml"
python "$SCRIPT_DIR/load_params.py" "$new_ncp_params_path"
export NCP_PARAMS_YML=$new_ncp_params_path

log info "Running all NCPs: Scenario ID=$NCP_RUN_SCENARIO_ID, Year=$NCP_RUN_YEAR"

## CAR - Carbon sequestration (R -> Python -> R)
CAR_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/CAR/tot_c_cur_$NCP_RUN_YEAR.tif"
if [ ! -f "$CAR_output_file" ]; then
    log info "Running CAR - Carbon sequestration"
    run_scripts_with_env "r" "$SCRIPT_DIR/CAR/1_CAR_S_CH.R"
    run_scripts_with_env "python" "$SCRIPT_DIR/CAR/2_CAR_S_CH.py"
    run_scripts_with_env "r" "$SCRIPT_DIR/CAR/3_CAR_S_CH.R"
else
    log info "Skipping CAR - Carbon sequestration as $CAR_output_file already exists"
fi

## FF - Food and feed (R only)
FF_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/FF/FF_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$FF_output_file" ]; then
    log info "Running FF - Crop production potential"
    run_scripts_with_env "r" "$SCRIPT_DIR/FF/1_FF_aggregate_masking.R"
else
    log info "Skipping FF - Crop production potential as $FF_output_file already exists"
fi

## HAB - Habitat quality (R -> Python)
HAB_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/HAB/$NCP_RUN_YEAR"
if [ ! -f "$HAB_output_folder/deg_sum_c.tif" ] || [ ! -f "$HAB_output_folder/quality_c.tif" ]; then
    log info "Running HAB - Habitat quality index"
    run_scripts_with_env "r" "$SCRIPT_DIR/HAB/0_HAB_threat_layers_generation.R"
    run_scripts_with_env "python" "$SCRIPT_DIR/HAB/1_HAB_S_CH.py"
else
    log info "Skipping HAB - Habitat quality index"
fi

## NDR - Nutrient Delivery Ratio (Python only)
NDR_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/NDR/$NCP_RUN_YEAR"
if [ ! -f "$NDR_output_folder/n_subsurface_export.tif" ] ||
   [ ! -f "$NDR_output_folder/n_surface_export.tif" ] ||
   [ ! -f "$NDR_output_folder/n_total_export.tif" ] ||
   [ ! -f "$NDR_output_folder/p_surface_export.tif" ]; then
    log info "Running NDR - Annual nutrient retention by vegetation"
    run_scripts_with_env "python" "$SCRIPT_DIR/NDR/1_NDR_S_CH.py"
else
    log info "Skipping NDR - Annual nutrient retention by vegetation"
fi

## POL - Pollination (Python -> R)
POL_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/POL/POL_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$POL_output_file" ]; then
    log info "Running POL - Habitat abundance for pollinators"
    run_scripts_with_env "python" "$SCRIPT_DIR/POL/1_POL_S_CH.py"
    run_scripts_with_env "r" "$SCRIPT_DIR/POL/2_POL_S_CH_aggregating.R"
else
    log info "Skipping POL - Habitat abundance for pollinators"
fi

## SDR - Sediment Delivery Ratio (Python only)
SDR_output_files="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/SDR/$NCP_RUN_YEAR/"
if [ ! -f "$SDR_output_files/avoided_erosion.tif" ] ||
   [ ! -f "$SDR_output_files/avoided_export.tif" ] ||
   [ ! -f "$SDR_output_files/sed_deposition.tif" ] ||
   [ ! -f "$SDR_output_files/sed_export.tif" ] ||
   [ ! -f "$SDR_output_files/stream.tif" ] ||
   [ ! -f "$SDR_output_files/rkls.tif" ] ||
   [ ! -f "$SDR_output_files/usle.tif" ]; then
    log info "Running SDR - Erosion control by sediment retention"
    run_scripts_with_env "python" "$SCRIPT_DIR/SDR/1_SDR_S_CH.py"
else
    log info "Skipping SDR - Erosion control by sediment retention"
fi

## WY - Water Yield (Python only)
WY_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/WY/$NCP_RUN_YEAR"
if [ ! -f "$WY_output_folder/watershed_results_wyield.shp" ] ||
   [ ! -f "$WY_output_folder/per_pixel/wyield.tif" ]; then
    log info "Running WY - Annual water yield"
    run_scripts_with_env "python" "$SCRIPT_DIR/WY/1_WY_S_CH.py"
else
    log info "Skipping WY - Annual water yield"
fi

## REC - Recreation potential (R only)
REC_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/REC/REC_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$REC_output_file" ]; then
    log info "Running REC - Recreation potential"
    run_scripts_with_env "r" "$SCRIPT_DIR/REC/1_REC.R"
else
    log info "Skipping REC - Recreation potential"
fi

log info "All NCPs have been run"
NCP_PARAMS_YML=$NCP_PARAMS_YML_BASE