#!/bin/bash

# Array of input files
input_files=(
    # 20 40
    "P_20_40_1.xml"
    # "P_20_40_2.xml"
    # "P_20_40_3.xml"
    # "P_20_40_4.xml"
    "P_20_40_old.xml"
    # 40 80
    "P_40_80_1.xml"
    # "P_40_80_2.xml"
    # "P_40_80_3.xml"
    # "P_40_80_4.xml"
    "P_40_80_old.xml"
    # 60 160
    "P_60_120_1.xml"
    # "P_60_120_2.xml"
    # "P_60_120_3.xml"
    # "P_60_120_4.xml"
    "P_60_120_old.xml"
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