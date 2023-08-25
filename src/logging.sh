#!/bin/bash

# Define a log file
LOG_FILE="my_log.txt"

# Define ANSI color codes
COLOR_RESET='\033[0m'      # Reset
COLOR_INFO='\033[0;32m'    # Green
COLOR_WARNING='\033[0;33m' # Yellow
COLOR_ERROR='\033[0;31m'   # Red

# Define host name
HOST=$(hostname)

# Custom logging function
# log <level> <message>
log() {
  local timestamp
  timestamp=$(date +"%Y-%m-%d %H:%M:%S")
  local log_message="[$timestamp] [$HOST] $2"

  # Check if a log level is provided (e.g., info, warning, error)
  if [[ ! -z "$1" ]]; then
    case "$1" in
      "info") log_message="${COLOR_INFO}$log_message${COLOR_RESET}" ;;
      "warning") log_message="${COLOR_WARNING}$log_message${COLOR_RESET}" ;;
      *) log_message="${COLOR_ERROR}$log_message${COLOR_RESET}" ;;
    esac
  fi

  echo -e "$log_message"
}