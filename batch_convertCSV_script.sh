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
MovingObstacle_new_summation_Re50_400_200_1rf_1prc
MovingObstacle_new_summation_Re100_400_200_1rf_1prc
MovingObstacle_new_summation_Re150_400_200_1rf_1prc
) 

# Loop through the hardcoded list of files
for filename in "${files[@]}"; do
    # Call the Python script for each file pattern
    echo "Processing files for pattern: $filename"
    pvpython convertCSV.py "$filename" $start $end
done