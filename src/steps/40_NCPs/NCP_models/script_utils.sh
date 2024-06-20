#!/bin/bash
# Additional functions for the NCP run scripts


: "
run_scripts - Run list of scripts sequentially

Known script types:
    - .R  - R
    - .py - Python

Parameters:
    $1 - List of paths to scripts

Returns:
    None

Raises:
    Logs an error if script type is not recognized
"
function run_scripts() {
    local scripts=("$@")
    for script in "${scripts[@]}"; do
        case "$script" in
            *.R)
                log debug "Running R script: $script"
#                Rscript "$script" # There were 16 warnings (use warnings() to see them)
                # to show warnings, use the following command
                Rscript --vanilla "$script"
                ;;
            *.py)
                log debug "Running Python script: $script"
                python "$script"
                ;;
            *)
                log error "Unknown script type: $script"
                ;;
        esac
    done
}
