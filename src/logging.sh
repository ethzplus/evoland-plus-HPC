#!/bin/bash

# Define ANSI color codes
COLOR_RESET='\033[0m'      # Reset
COLOR_DEBUG='\033[0;34m'   # Blue
COLOR_INFO='\033[0;32m'    # Green
COLOR_WARNING='\033[0;33m' # Yellow
COLOR_ERROR='\033[0;31m'   # Red

# Define host name
HOST=$(hostname)
export FUTURE_EI_LOG_LEVEL=${FUTURE_EI_LOG_LEVEL:-info} # default log level is info

# get numeric representation of log level (debug: 0, info: 1, warning: 2, error: 3)
level_to_number() {
  case "$1" in
    "debug") echo 0 ;;
    "info") echo 1 ;;
    "warning") echo 2 ;;
    *) echo 3 ;;
  esac
}

FUTURE_EI_LOG_LEVEL_NUMBER=$(level_to_number "$FUTURE_EI_LOG_LEVEL")

# Custom logging function
# log <level> <message>
log() {
  local timestamp
  timestamp=$(date +"%Y-%m-%d %H:%M:%S")
  local log_message="[$timestamp] [$HOST] $2"
  # message level is lower than log level, leave function
  if [[ $(level_to_number "$1") -lt $FUTURE_EI_LOG_LEVEL_NUMBER ]]; then
    return
  fi

  # Check if a log level is provided (e.g., debug, info, warning, error)
  if [[ -n "$1" ]]; then
    case "$1" in
      "debug") log_message="${COLOR_DEBUG}$log_message${COLOR_RESET}" ;;
      "info") log_message="${COLOR_INFO}$log_message${COLOR_RESET}" ;;
      "warning") log_message="${COLOR_WARNING}$log_message${COLOR_RESET}" ;;
      *) log_message="${COLOR_ERROR}$log_message${COLOR_RESET}" ;;
    esac
  fi

  echo -e "$log_message"
}