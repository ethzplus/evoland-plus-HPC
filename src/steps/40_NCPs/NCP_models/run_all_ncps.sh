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
CAR_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/CAR/tot_c_cur_$NCP_RUN_YEAR.tif"
if [ ! -f "$CAR_output_file" ]; then
    log info "Running CAR - Carbon sequestration"
    run_scripts "$SCRIPT_DIR/CAR/1_CAR_S_CH.R" \
                "$SCRIPT_DIR/CAR/2_CAR_S_CH.py" \
                "$SCRIPT_DIR/CAR/3_CAR_S_CH.R"
else
    log info "Skipping CAR - Carbon sequestration as $CAR_output_file already exists"
fi

## FF - Food and feed
## Indicator: Crop production potential (`ecocrop`)
FF_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/FF/FF_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$FF_output_file" ]; then
    log info "Running FF - Crop production potential"
    # needs FF/data_preparation/0_monthly_average_dataprep.R from prepare_ncps.sh
    run_scripts "$SCRIPT_DIR/FF/1_FF_aggregate_masking.R"
else
    log info "Skipping FF - Crop production potential as $FF_output_file already exists"
fi

## HAB - Habitat creation and maintenance
# check natcap version for quickfix hint if needed
# (natcap.invest.__version__ < 3.14.0)
show_hint=$(python -c "from natcap.invest import __version__; from packaging.version import Version; print(Version(__version__) < Version('3.14.0'))")
if [ "$show_hint" = "True" ]; then
    log warning "As you have a natcap version < 3.14.0, you might need to apply a quickfix to the package source code.
    Find it in 'src/steps/40_NCPs/NCP_models/HAB/1_HAB_S_CH.py'."
fi
## Indicator: Habitat quality index
HAB_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/HAB/$NCP_RUN_YEAR"
if [ ! -f "$HAB_output_folder/deg_sum_c.tif" ] || [ ! -f "$HAB_output_folder/quality_c.tif" ]; then
    log info "Running HAB - Habitat quality index"
    run_scripts "$SCRIPT_DIR/HAB/0_HAB_threat_layers_generation.R" \
                "$SCRIPT_DIR/HAB/1_HAB_S_CH.py"
else
    log info "Skipping HAB - Habitat quality index as $HAB_output_folder/deg_sum_c.tif and $HAB_output_folder/quality_c.tif already exist"
fi

## NDR - Nutrient Delivery Ratio
## Indicator: Annual nutrient retention by vegetation
NDR_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/NDR/$NCP_RUN_YEAR"
if [ ! -f "$NDR_output_folder/n_subsurface_export.tif" ] ||
   [ ! -f "$NDR_output_folder/n_surface_export.tif" ] ||
   [ ! -f "$NDR_output_folder/n_total_export.tif" ] ||
   [ ! -f "$NDR_output_folder/p_surface_export.tif" ]; then
    log info "Running NDR - Annual nutrient retention by vegetation"
    run_scripts "$SCRIPT_DIR/NDR/1_NDR_S_CH.py"
else
    log info "Skipping NDR - Annual nutrient retention by vegetation as $NDR_output_folder/n_subsurface_export.tif and $NDR_output_folder/n_surface_export.tif and $NDR_output_folder/n_total_export.tif and $NDR_output_folder/p_surface_export.tif already exist"
fi

## POL - Pollination and dispersal of seeds
## Indicator: Habitat abundance for pollinators
POL_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/POL/POL_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$POL_output_file" ]; then
    log info "Running POL - Habitat abundance for pollinators"
    run_scripts "$SCRIPT_DIR/POL/1_POL_S_CH.py" \
                "$SCRIPT_DIR/POL/2_POL_S_CH_aggregating.R"
else
    log info "Skipping POL - Habitat abundance for pollinators as $POL_output_file already exists"
fi

## SDR - Formation, protection and decontamination of soils
## Indicator: Erosion control by sediment retention
SDR_output_files="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/SDR/$NCP_RUN_YEAR/"
if [ ! -f "$SDR_output_files/avoided_erosion.tif" ] ||
   [ ! -f "$SDR_output_files/avoided_export.tif" ] ||
   [ ! -f "$SDR_output_files/sed_deposition.tif" ] ||
   [ ! -f "$SDR_output_files/sed_export.tif" ] ||
   [ ! -f "$SDR_output_files/stream.tif" ] ||
   [ ! -f "$SDR_output_files/rkls.tif" ]
   [ ! -f "$SDR_output_files/usle.tif" ]; then
    log info "Running SDR - Erosion control by sediment retention"
    run_scripts "$SCRIPT_DIR/SDR/1_SDR_S_CH.py"
else
    log info "Skipping SDR - Erosion control by sediment retention as
    $SDR_output_files/avoided_erosion.tif, avoided_export.tif, sed_deposition.tif, sed_export.tif, stream.tif, rkls.tif and usle.tif already exist"
fi

## WY - Regulation of freshwater quantity, location and timing
## Indicator: Annual water yield
WY_output_folder="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/WY/$NCP_RUN_YEAR"
if [ ! -f "$WY_output_folder/watershed_results_wyield.shp" ] ||
   [ ! -f "$WY_output_folder/watershed_results_wyield.shx" ] ||
   [ ! -f "$WY_output_folder/watershed_results_wyield.dbf" ] ||
   [ ! -f "$WY_output_folder/watershed_results_wyield.prj" ] ||
   [ ! -f "$WY_output_folder/watershed_results_wyield.csv" ] ||
   [ ! -f "$WY_output_folder/subwatershed_results_wyield.shp" ] ||
   [ ! -f "$WY_output_folder/subwatershed_results_wyield.shx" ] ||
   [ ! -f "$WY_output_folder/subwatershed_results_wyield.dbf" ] ||
   [ ! -f "$WY_output_folder/subwatershed_results_wyield.prj" ] ||
   [ ! -f "$WY_output_folder/subwatershed_results_wyield.csv" ] ||
   [ ! -f "$WY_output_folder/per_pixel/aet.tif" ] ||
   [ ! -f "$WY_output_folder/per_pixel/fractp.tif" ] ||
   [ ! -f "$WY_output_folder/per_pixel/wyield.tif" ]; then
    log info "Running WY - Annual water yield"
    run_scripts "$SCRIPT_DIR/WY/1_WY_S_CH.py"
else
    log info "Skipping WY - Annual water yield as $WY_output_folder/watershed_results_wyield.* and $WY_output_folder/subwatershed_results_wyield.* and $WY_output_folder/per_pixel already exist"
fi

## REC - Recreation potential
## Indicator: Recreation potential (RP) provided by ecosystems
REC_output_file="$NCP_RUN_OUTPUT_DIR/$NCP_RUN_SCENARIO_ID/REC/REC_S_CH_$NCP_RUN_YEAR.tif"
if [ ! -f "$REC_output_file" ]; then
    log info "Running REC - Recreation potential"
    run_scripts "$SCRIPT_DIR/REC/1_REC.R"
else
    log info "Skipping REC - Recreation potential as $REC_output_file already exists"
fi

# ------------------------------------------------------------------------------

log info "All NCPs have been run"

NCP_PARAMS_YML=$NCP_PARAMS_YML_BASE