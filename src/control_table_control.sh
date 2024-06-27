#!/bin/bash
# ----------------------------------------------
# Distribute list of rows to SLURM jobs
# Functions to split and unify the control table
# ----------------------------------------------
: "
get_rows - Function to distribute list of rows to SLURM jobs

Parameters:
    $1 - Current task ID
    $2 - Total number of tasks/array jobs
    $3 - Minimum task ID
    $4 - Total number of scenarios

Returns:
    A string with the start and end row numbers separated by a space.

Raises:
    Logs an error and returns if the number of tasks/array jobs is larger than the number of scenarios.
"
function get_rows() {
    # num_sim needs to be >= task_count
    if [ "$4" -lt "$2" ]; then
        log error "Number of tasks/array jobs ($2) cannot be larger than the number of scenarios ($4)."
        return
    fi
    # Calculate the number of full rows each task gets
    full_rows=$(( $4 / $2 ))
    # Calculate the remainder of the division
    remainder=$(( $4 % $2 ))
    # Calculate start_row and end_row
    start_row=$(( (($1 - $3) * full_rows) + 2 + ($1 - $3 <= remainder ? ($1 - $3) : remainder) ))
    end_row=$(( start_row + full_rows + ($1 - $3 < remainder ? 1 : 0) ))
    echo "$start_row $end_row"
}

: "
split_control_table - Function to split the control table

Save a chunk of the control table to a temporary directory.

Uses the environment variables:
    LULCC_M_SIM_CONTROL_TABLE      - Name of the simulation control table
    LULCC_M_SIM_CONTROL_TABLE_FULL - Name of the full simulation control table
    SLURM_ARRAY_JOB_ID             - Job array ID number
    SLURM_ARRAY_TASK_ID            - Job array ID (index) number
    SLURM_ARRAY_TASK_COUNT         - Total number of tasks in a job array
    SLURM_ARRAY_TASK_MIN           - Job array's minimum ID (index) number
    LULCC_START_ROW                - Start row number
    LULCC_END_ROW                  - End row number
"
function split_control_table() {
    # Show used LULCC_M_SIM_CONTROL_TABLE file
    log info "Splitting control table at $LULCC_M_SIM_CONTROL_TABLE"
    # number of rows in simulation control file (minus header, ignoring empty
    # lines in the end)
    num_sims=$(($(grep -c -v '^$' "$LULCC_M_SIM_CONTROL_TABLE")-1))

    # use function get_rows() to get start_row and end_row
    ## if SBATCH --array=1-4 => SLURM_ARRAY_TASK_ID=1 through 4;
    # SLURM_ARRAY_TASK_COUNT=4; SLURM_ARRAY_TASK_MIN=1
    log debug "task_id: $SLURM_ARRAY_TASK_ID, task_count: $SLURM_ARRAY_TASK_COUNT, task_min: $SLURM_ARRAY_TASK_MIN"
    start_end_rows=$(get_rows "$SLURM_ARRAY_TASK_ID" "$SLURM_ARRAY_TASK_COUNT" "$SLURM_ARRAY_TASK_MIN" "$num_sims")
    log debug "num_sims: $num_sims, start_end_rows: $start_end_rows"

    LULCC_START_ROW=$(echo "$start_end_rows" | cut -d' ' -f1)
    LULCC_END_ROW=$(echo "$start_end_rows" | cut -d' ' -f2)
    export LULCC_START_ROW
    export LULCC_END_ROW

    # save subset of simulation control file to temporary directory
    log debug "Saving subset of simulation control file to temporary directory"
    new_control_table_path="$LULCC_CH_HPC_DIR/Tools/Simulation_control_subsets/Simulation_control-${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID.csv"
    mkdir -p "$(dirname "$new_control_table_path")"
    head -n 1 "$LULCC_M_SIM_CONTROL_TABLE" > "$new_control_table_path"
    sed -n "$LULCC_START_ROW,$((LULCC_END_ROW-1))"p "$LULCC_M_SIM_CONTROL_TABLE" >> "$new_control_table_path"
    # remember the path to the full simulation control file
    export LULCC_M_SIM_CONTROL_TABLE_FULL
    LULCC_M_SIM_CONTROL_TABLE_FULL="$LULCC_M_SIM_CONTROL_TABLE"
    # set the path to the subset of the simulation control file
    LULCC_M_SIM_CONTROL_TABLE="$new_control_table_path"

    log info "Running rows [$LULCC_START_ROW, ${LULCC_END_ROW}[ of $num_sims"
}

: "
merge_control_table - Function to merge the control table

Merge the temporary control table chunks back into the original control table.

Uses the environment variables:
    LULCC_M_SIM_CONTROL_TABLE      - Name of the simulation control table
    LULCC_M_SIM_CONTROL_TABLE_FULL - Name of the full simulation control table
    TMPDIR                         - Path to the temporary directory
    LULCC_START_ROW                - Start row number
    LULCC_END_ROW                  - End row number
"
function merge_control_table() {
    # Show used LULCC_M_SIM_CONTROL_TABLE file
    log info "Merging control table at $LULCC_M_SIM_CONTROL_TABLE into $LULCC_M_SIM_CONTROL_TABLE_FULL"
    # Use lockfile to prevent multiple jobs from writing to the same file
    lockfile="$LULCC_M_SIM_CONTROL_TABLE_FULL.lock"
    # Wait for lockfile to be available
    while [ -f "$lockfile" ]; do
        log debug "Waiting for lockfile $lockfile to be available"
        sleep 5
    done
    # Create lockfile
    touch "$lockfile"
    log debug "Created lockfile $lockfile"
    # Set a trap to remove the lockfile on exit
    trap 'rm -f "$lockfile"; echo "Trap: Removed lockfile $lockfile"; exit $?' INT TERM EXIT

    # Create new output file if it does not exist, so processes can work in
    # parallel and concurrently
    # at $FUTURE_EI_OUTPUT_DIR/basename($LULCC_M_SIM_CONTROL_TABLE_FULL)
    dest_table="$FUTURE_EI_OUTPUT_DIR/$(basename "$LULCC_M_SIM_CONTROL_TABLE_FULL")"
    if [ ! -f "$dest_table" ]; then
        cp "$LULCC_M_SIM_CONTROL_TABLE_FULL" "$dest_table"
    fi
    LULCC_M_SIM_CONTROL_TABLE_FULL="$dest_table"
    # Replace the extracted rows in the full simulation control file with the
    # rows from the temporary simulation control file
    # 1. Delete the rows from LULCC_START_ROW to LULCC_END_ROW-1
    sed -i "$LULCC_START_ROW,$((LULCC_END_ROW-1))"d "$LULCC_M_SIM_CONTROL_TABLE_FULL"
    # sed -i: edit file in place
    #      d: delete the specified lines (from start_row to end_row-1)
    # 2. delete header in temporary simulation control file
    sed -i 1d "$LULCC_M_SIM_CONTROL_TABLE"
    # sed -i: edit file in place
    #     1d: delete the first line (header)
    # 3. copy all rows from LULCC_M_SIM_CONTROL_TABLE to
    # LULCC_M_SIM_CONTROL_TABLE_FULL at the corresponding position
    sed -i "$((LULCC_START_ROW-1)) r $LULCC_M_SIM_CONTROL_TABLE" "$LULCC_M_SIM_CONTROL_TABLE_FULL"
    # sed -i: edit file in place
    #      r: read file and append to the line (start_row-1)

    # Remove lockfile
    rm "$lockfile"
    # show folder where lockfile was created
    if [ -f "$lockfile" ]; then
        log error "Could not remove lockfile $lockfile"
    fi

    # Remove temporary simulation control file
    rm "$LULCC_M_SIM_CONTROL_TABLE"
    LULCC_M_SIM_CONTROL_TABLE="$LULCC_M_SIM_CONTROL_TABLE_FULL"

    log info "Merged rows [$LULCC_START_ROW, ${LULCC_END_ROW}[ of $num_sims"
}
export -f split_control_table
export -f merge_control_table
