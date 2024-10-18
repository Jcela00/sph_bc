#!/bin/bash

# Check if exactly 2 arguments (start and end) are provided
if [ "$#" -ne 2 ]; then
    echo "You must enter exactly 2 command line arguments: start and end"
    exit 1
fi

# Get the start and end arguments
start=$1
end=$2

# List of filename beginings
files=(
    "TaylorCouette_new_summation_10_10_1rf_1prc"
    "TaylorCouette_new_summation_10_10_2rf_1prc"
    "TaylorCouette_new_summation_10_10_3rf_1prc"
    "TaylorCouette_old_summation_10_10_1rf_1prc"
    "TaylorCouette_new_summation_15_15_1rf_1prc"
    "TaylorCouette_new_summation_15_15_2rf_1prc"
    "TaylorCouette_new_summation_15_15_3rf_1prc"
    "TaylorCouette_old_summation_15_15_1rf_1prc"
    "TaylorCouette_new_summation_20_20_1rf_1prc"
    "TaylorCouette_new_summation_20_20_2rf_1prc"
    "TaylorCouette_new_summation_20_20_3rf_1prc"
    "TaylorCouette_old_summation_20_20_1rf_1prc"
    "TaylorCouette_new_summation_30_30_1rf_1prc"
    "TaylorCouette_new_summation_30_30_2rf_1prc"
    "TaylorCouette_new_summation_30_30_3rf_1prc"
)

# Loop through the hardcoded list of files
for filename in "${files[@]}"; do
    # Call the Python script for each file pattern
    echo "Processing files for pattern: $filename"
    pvpython convertCSV.py "$filename" $start $end
done