#!/bin/bash

# Array of input files
input_files=(
XML/Flower/Flower_old_60.xml
XML/Flower/Flower_new_60.xml
XML/Flower/Flower_old_120.xml
XML/Flower/Flower_new_120.xml
XML/Flower/Flower_old_240.xml
XML/Flower/Flower_new_240.xml
)

# Path to the simulation executable
executable="./sph_dlb"

# Loop over the input files and run the simulation sequentially
for input_file in "${input_files[@]}"; do
    echo "Running simulation with input file: $input_file"
    $executable $input_file

    echo "Simulation for $input_file completed."
done

echo "All simulations completed."