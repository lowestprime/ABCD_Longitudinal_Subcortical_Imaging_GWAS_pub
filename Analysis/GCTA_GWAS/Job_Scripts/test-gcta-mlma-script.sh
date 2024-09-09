#!/bin/bash

# Test script for GCTA MLMA parallelized script

set -euo pipefail

# Source the main script to access its functions
source ./improved-gcta-mlma-script.sh

# Test 1: Check if required software is available
test_software_availability() {
    echo "Testing software availability..."
    if command -v parallel &> /dev/null && [[ -x "$gcta" ]]; then
        echo "PASS: Required software is available"
    else
        echo "FAIL: Required software is not available"
        return 1
    fi
}

# Test 2: Check if directories exist
test_directory_structure() {
    echo "Testing directory structure..."
    local dirs=(
        "$base_dir"
        "$pheno_dir"
        "$covar_dir"
        "$results_dir"
    )
    for dir in "${dirs[@]}"; do
        if [[ ! -d "$dir" ]]; then
            echo "FAIL: Directory $dir does not exist"
            return 1
        fi
    done
    echo "PASS: All required directories exist"
}

# Test 3: Check file existence for a sample combination
test_file_existence() {
    echo "Testing file existence..."
    local test_pop="EUR"
    local test_sex="F"
    local test_pheno="smri_vol_scs_caudalanteriorcingulate_ROC0_2"
    if check_files_exist "$test_pop" "$test_sex" "$test_pheno"; then
        echo "PASS: All required files exist for $test_pop $test_sex $test_pheno"
    else
        echo "FAIL: Missing files for $test_pop $test_sex $test_pheno"
        return 1
    fi
}

# Test 4: Check skip combination logic
test_skip_combination() {
    echo "Testing skip combination logic..."
    if should_skip "EUR" "M" "smri_vol_scs_wholeb_ROC0_2"; then
        echo "PASS: Correctly identified combination to skip"
    else
        echo "FAIL: Failed to identify combination to skip"
        return 1
    fi
    if ! should_skip "EUR" "F" "smri_vol_scs_caudalanteriorcingulate_ROC0_2"; then
        echo "PASS: Correctly identified combination not to skip"
    else
        echo "FAIL: Incorrectly identified combination to skip"
        return 1
    fi
}

# Test 5: Run a small-scale GCTA MLMA test
test_gcta_mlma() {
    echo "Testing GCTA MLMA execution..."
    local test_pop="EUR"
    local test_sex="F"
    local test_pheno="smri_vol_scs_caudalanteriorcingulate_ROC0_2"
    local test_scratch_dir="/tmp/gcta_test_${RANDOM}"
    mkdir -p "$test_scratch_dir"
    if run_gcta_mlma "$test_pop" "$test_sex" "$test_pheno" "$test_scratch_dir"; then
        echo "PASS: GCTA MLMA test run completed successfully"
        rm -rf "$test_scratch_dir"
    else
        echo "FAIL: GCTA MLMA test run failed"
        rm -rf "$test_scratch_dir"
        return 1
    fi
}

# Run all tests
run_all_tests() {
    local tests=(
        test_software_availability
        test_directory_structure
        test_file_existence
        test_skip_combination
        test_gcta_mlma
    )
    local failed=0
    for test in "${tests[@]}"; do
        if ! $test; then
            ((failed++))
        fi
        echo
    done
    if ((failed == 0)); then
        echo "All tests passed successfully!"
    else
        echo "$failed test(s) failed."
        return 1
    fi
}

# Execute all tests
run_all_tests