#!/bin/bash
# --------------------------------------------------------------------------
# This script is used for cleaning up and saving the output of the future EI
# runs. It checks if the output folder exists and creates it if not.
# It backs up any previous output folder and moves the output folders to the
# new output folder.
# --------------------------------------------------------------------------

# List of folders to move (variable name)
DIRS_TO_MOVE=(
  "LULCC_CH_OUTPUT_BASE_DIR"
  "NSDM_OUTPUT_BASE_DIR"
  "SPAGG_OUTPUT_BASE_DIR"
  "NCP_OUTPUT_BASE_DIR"
)

: "
backup_dir - Function to backup a directory

This function takes a directory as an argument, backs it up by moving it to a
new location $dir-backup and then creates a new directory at the original
location. If the backup directory already exists, it is deleted.
Arguments:
  dir: The directory to backup
Returns:
  None
"
backup_dir() {
  local dir="$1"
  local backup_dir="${dir}-backup"

  log info "Backing up $dir to $backup_dir"

  # If backup folder exists, delete it
  if [ -d "$backup_dir" ]; then
    log debug "Backup folder already exists, deleting it"
    rm -r "$backup_dir"
  fi

  # Move dir to backup_dir
  mv "$dir" "$backup_dir"
}

: "
move_dir - Function to move a directory to a new location

This function takes a source directory and a destination directory as
arguments. It checks if the source directory is set and if it exists,
then moves it to the destination directory.

Arguments:
  src: The source directory to move, unresolved variable
  dest: The destination directory to move the source directory to
Returns:
  None
"
move_dir() {
  local src="$1"
  local dest="$2"

  # check if variable is set
  if [ -z "${!src}" ]; then
    log warning "Folder $src is not set"
    return
  fi

  # check if folder exists
  if [ ! -d "${!src}" ]; then
    log warning "Folder ${!src} does not exist"
    return
  fi

  # move folder to $dest
  log info "Moving $src (${!src}) to $dest"
  mv "${!src}" "$dest"
}

# Check if $FUTURE_EI_OUTPUT_DIR exists and create it if not
# Backup old folder if it exists
if [ -d "$FUTURE_EI_OUTPUT_DIR" ]; then
  backup_dir "$FUTURE_EI_OUTPUT_DIR"
fi

# Create new output folder
log info "Creating $FUTURE_EI_OUTPUT_DIR"
mkdir -p "$FUTURE_EI_OUTPUT_DIR" # -p creates parent directories if they don't exist

# Move folders - output folder is created and thus empty
for var in "${DIRS_TO_MOVE[@]}"; do
  move_dir "$var" "$FUTURE_EI_OUTPUT_DIR"
done