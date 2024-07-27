#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments"
    exit 1
fi

# Get the arguments
filename=$1
start=$2
end=$3

# Execute the command
pvpython convertCSV.py $filename $start $end