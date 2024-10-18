#!/bin/bash

# Array of input files
input_files=(
    "TC_10_1.xml"   
    "TC_10_2.xml"
    "TC_10_3.xml"
    "TC_10_old.xml"   
    "TC_15_1.xml"   
    "TC_15_2.xml"
    "TC_15_3.xml"
    "TC_15_old.xml"   
    "TC_20_1.xml"   
    "TC_20_2.xml"
    "TC_20_3.xml"
    "TC_20_old.xml"
    "TC_30_1.xml"
    "TC_30_2.xml"
    "TC_30_3.xml"
    "TC_30_old.xml"
    "TC_40_1.xml"
    "TC_40_2.xml"
    "TC_40_3.xml"
    "TC_40_old.xml"
    # Add more file names here
)

# Number of processes to use
num_procs=1

# Path to the simulation executable
executable="./sph_dlb"

# Loop over the input files and run the simulation sequentially
for input_file in "${input_files[@]}"; do
    echo "Running simulation with input file: $input_file"
    mpirun -np $num_procs $executable XMLs/XML_TC/$input_file

    echo "Simulation for $input_file completed."
done

echo "All simulations completed."