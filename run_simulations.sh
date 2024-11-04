#!/bin/bash

# Array of input files
input_files=(
  # Poiseuille_old_30.xml
  # Poiseuille_old_60.xml
  # Poiseuille_old_90.xml
  # Poiseuille_old_120.xml
  # Poiseuille_old_180.xml
  # Poiseuille_old_240.xml

  # Poiseuille_30.xml
  # Poiseuille_60.xml
  # Poiseuille_90.xml
  # Poiseuille_120.xml
  Poiseuille_180.xml
  Poiseuille_240.xml

  # TC_10_1.xml
  # TC_20_1.xml
  # TC_30_1.xml
  # TC_40_1.xml
  # TC_50_1.xml

  # TC_10_2.xml
  # TC_20_2.xml
  # # TC_30_2.xml
  # # TC_40_2.xml
  # TC_50_2.xml

  # # TC_10_1_sum.xml
  # # TC_20_1_sum.xml
  # TC_30_1_sum.xml
  # TC_40_1_sum.xml
  # TC_50_1_sum.xml

  # # TC_old_10.xml
  # # TC_old_20.xml
  # TC_old_30.xml
  # TC_old_40.xml
  # TC_old_50.xml

  # TC_10_0.xml
  # TC_20_0.xml
  # TC_30_0.xml
  # TC_40_0.xml
  # TC_50_0.xml

)

# Number of processes to use
num_procs=1

# Path to the simulation executable
executable="./sph_dlb"

# Loop over the input files and run the simulation sequentially
for input_file in "${input_files[@]}"; do
    echo "Running simulation with input file: $input_file"
    $executable XML/poiseuilleDiff/$input_file

    echo "Simulation for $input_file completed."
done

echo "All simulations completed."