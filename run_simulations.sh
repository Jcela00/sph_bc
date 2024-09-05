#!/bin/bash

# Array of input files
input_files=(
    "TC_40_40_1.xml"   
    "TC_40_40_2.xml"
    "TC_40_40_3.xml"
    "TC_40_40_old.xml"   
    "TC_80_80_1.xml"   
    "TC_80_80_2.xml"
    "TC_80_80_3.xml"
    "TC_80_80_old.xml"
    "TC_120_120_1.xml"
    "TC_120_120_2.xml"
    "TC_120_120_3.xml"
    "TC_120_120_old.xml"
    "TC_160_160_1.xml"
    "TC_160_160_2.xml"
    "TC_160_160_3.xml"
    "TC_160_160_old.xml"
    "TC_240_240_1.xml"
    "TC_240_240_2.xml"
    "TC_240_240_3.xml"
    "TC_240_240_old.xml"
    # Add more file names here
)

# Number of processes to use
num_procs=4

# Path to the simulation executable
executable="./sph_dlb"

# Loop over the input files and run the simulation sequentially
for input_file in "${input_files[@]}"; do
    echo "Running simulation with input file: $input_file"
    mpirun -np $num_procs $executable $input_file

    echo "Simulation for $input_file completed."
done

echo "All simulations completed."