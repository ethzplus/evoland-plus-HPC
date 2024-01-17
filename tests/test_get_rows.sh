#!/bin/bash
# -------------------------------
# Tests for the get_rows function
# -------------------------------
source ../src/bash_common.sh
source ../src/get_rows.sh

# Variable to keep track of whether any tests have failed
tests_failed=0

: "
run_test - Function to run a test and compare the expected output with the actual output

Parameters:
    $1 - Expected output
    $2 - Actual output
    $3 - Optional: Error expected (0: no, 1: yes)

Returns:
    The expected output if the test passes, otherwise nothing.
"
function run_test() {
    local expected_output=$1
    local output=$2
    local error_expected=${3:-0}

    # Check if an error was expected
    if [ "$error_expected" -eq 1 ]; then
      # Output does not contain expected error message
      if ! grep -q "Number of tasks/array jobs" <<< "$output"; then
          echo "Test failed: expected error but got $output"
          tests_failed=1
      fi
    else
        # Output is not equal to expected output
        if [ "$output" != "$expected_output" ]; then
            echo "Test failed: expected $expected_output but got $output"
            tests_failed=1
        fi
    fi
}

# Test 1 - no array, 1 scenario
run_test "2 3" "$(get_rows 1 1 1 1)"
# Test 2 - no array, 2 scenarios
run_test "2 4" "$(get_rows 1 1 1 2)"
# Test 3 - no array, 100 scenarios
run_test "2 102" "$(get_rows 1 1 1 100)"

# Test 4 - 2 array jobs, 1 scenario - error: task_count > num_sims
run_test "" "$(get_rows 1 2 1 1)" 1

# Test 5 - 2 array jobs, 2 scenarios
run_test "2 3" "$(get_rows 1 2 1 2)"
run_test "3 4" "$(get_rows 2 2 1 2)"
# Test 6 - 2 array jobs, 3 scenarios
run_test "2 4" "$(get_rows 1 2 1 3)"
run_test "4 5" "$(get_rows 2 2 1 3)"
# Test 7 - 2 array jobs, 4 scenarios
run_test "2 4" "$(get_rows 1 2 1 4)"
run_test "4 6" "$(get_rows 2 2 1 4)"

# Test 7 - 2 array jobs, 4 scenarios, different start
run_test "2 4" "$(get_rows 0 2 0 4)"
run_test "4 6" "$(get_rows 1 2 0 4)"

# Test 8 - 3 array jobs, 2 scenario - error
run_test "" "$(get_rows 0 3 0 2)" 1

# Test 9 - 3 array jobs, 3 scenarios
run_test "2 3" "$(get_rows 1 3 1 3)"
run_test "3 4" "$(get_rows 2 3 1 3)"
run_test "4 5" "$(get_rows 3 3 1 3)"

# Test 10 - 3 array jobs, 30 scenarios
run_test "2 12" "$(get_rows 1 3 1 30)"
run_test "12 22" "$(get_rows 2 3 1 30)"
run_test "22 32" "$(get_rows 3 3 1 30)"
# Test 11 - 3 array jobs, 31 scenarios
run_test "2 13" "$(get_rows 1 3 1 31)"
run_test "13 23" "$(get_rows 2 3 1 31)"
run_test "23 33" "$(get_rows 3 3 1 31)"
# Test 12 - 3 array jobs, 32 scenarios
run_test "2 13" "$(get_rows 1 3 1 32)"
run_test "13 24" "$(get_rows 2 3 1 32)"
run_test "24 34" "$(get_rows 3 3 1 32)"

# Test 13 - 7 array jobs, 100 scenarios, different start
run_test "2 17" "$(get_rows 0 7 0 100)" # 15 scenarios
run_test "17 32" "$(get_rows 1 7 0 100)" # 15 scenarios
run_test "32 46" "$(get_rows 2 7 0 100)" # 14 scenarios
run_test "46 60" "$(get_rows 3 7 0 100)" # 14 scenarios
run_test "60 74" "$(get_rows 4 7 0 100)" # 14 scenarios
run_test "74 88" "$(get_rows 5 7 0 100)" # 14 scenarios
run_test "88 102" "$(get_rows 6 7 0 100)" # 14 scenarios
# Test 14 - 7 array jobs, 104 scenarios, different start
run_test "2 17" "$(get_rows 0 7 0 104)" # 15 scenarios
run_test "17 32" "$(get_rows 1 7 0 104)" # 15 scenarios
run_test "32 47" "$(get_rows 2 7 0 104)" # 15 scenarios
run_test "47 62" "$(get_rows 3 7 0 104)" # 15 scenarios
run_test "62 77" "$(get_rows 4 7 0 104)" # 15 scenarios
run_test "77 92" "$(get_rows 5 7 0 104)" # 15 scenarios
run_test "92 106" "$(get_rows 6 7 0 104)" # 14 scenarios

# Test 15 - 4 array jobs, 100 scenarios, different start
run_test "2 27" "$(get_rows -1 4 -1 100)"

# Test 16 - 4 array jobs, 100 scenarios, different start (unrealistic)
run_test "-24 2" "$(get_rows 0 4 1 100)"

# Test 17 - no array, 0 scenarios - error: task_count > num_sims
run_test "" "$(get_rows 1 1 1 0)" 1

# Test 18 - no array, 1e6 scenarios
run_test "2 1000002" "$(get_rows 1 1 1 1000000)"


# If any tests have failed, exit with a non-zero exit code
if [ "$tests_failed" -eq 1 ]; then
    log error "Some tests failed, see above for details."
    exit 1
fi