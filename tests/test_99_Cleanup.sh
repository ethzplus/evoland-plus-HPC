#!/bin/bash
# ------------------------------------------------
# Tests for the 99_Cleanup.sh functions and script
# ------------------------------------------------
source ../src/bash_common.sh

# Define the directories to be moved
DIRS_TO_MOVE=(
  "LULCC_CH_OUTPUT_BASE_DIR"
  "NSDM_OUTPUT_BASE_DIR"
  "SPAGG_OUTPUT_BASE_DIR"
  "NCP_OUTPUT_BASE_DIR"
)

# Use /tmp as the temporary directory
export tmp
tmp="/tmp"

# Define the output directory
FUTURE_EI_INPUT_DIR="$tmp/future-ei-input-test"
export FUTURE_EI_OUTPUT_DIR
FUTURE_EI_OUTPUT_DIR="$tmp/future-ei-output-test"
log debug "Overwriting FUTURE_EI_OUTPUT_DIR: $FUTURE_EI_OUTPUT_DIR"

# Create the variables and directories to be moved
for dir in "${DIRS_TO_MOVE[@]}"; do
  declare "$dir"="$FUTURE_EI_INPUT_DIR/$dir"
  log debug "Overwriting $dir: ${!dir}"
  mkdir -p "${!dir}"
  touch "${!dir}/test-file-in-$dir"
done

# Create existing output directory to be backed up
mkdir -p "$FUTURE_EI_OUTPUT_DIR"
touch "$FUTURE_EI_OUTPUT_DIR/to-be-backed-up"

# Run the cleanup script
source ../src/steps/99_Cleanup.sh

# Check if the directories were moved correctly
for dir in "${DIRS_TO_MOVE[@]}"; do
  moved_dir="$FUTURE_EI_OUTPUT_DIR/$dir"
  if [ ! -d "$moved_dir" ]; then
    echo "Test failed: $moved_dir does not exist"
    exit 1
  fi
  if [ ! -f "$moved_dir/test-file-in-$dir" ]; then
    echo "Test failed: $moved_dir/test-file-in-$dir does not exist"
    exit 1
  fi
done

# Check if the output directory was backed up
if [ ! -d "$FUTURE_EI_OUTPUT_DIR-backup" ]; then
  echo "Test failed: $FUTURE_EI_OUTPUT_DIR-backup does not exist"
  exit 1
fi
if [ ! -f "$FUTURE_EI_OUTPUT_DIR-backup/to-be-backed-up" ]; then
  echo "Test failed: $FUTURE_EI_OUTPUT_DIR-backup/to-be-backed-up does not exist"
  exit 1
fi

# Clean up
rm -r "$FUTURE_EI_OUTPUT_DIR"
rm -r "$FUTURE_EI_OUTPUT_DIR-backup"
rm -r "$FUTURE_EI_INPUT_DIR"

log info "All tests passed"
exit 0