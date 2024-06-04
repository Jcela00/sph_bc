#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "You must enter exactly 2 command line arguments"
    exit 1
fi

# Get the arguments
filename=$1
count=$2

# Execute the command
pvpython convertCSV.py $filename $count