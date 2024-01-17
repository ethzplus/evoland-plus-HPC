#!/bin/bash
# -------------------------------------------------
# Function to distribute list of rows to SLURM jobs
# -------------------------------------------------
: "
get_rows - Function to distribute list of rows to SLURM jobs

Parameters:
    $1 - Task ID
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

export -f get_rows