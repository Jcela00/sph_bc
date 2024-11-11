#!/bin/bash

# Array of input files
input_files=(

  # Poiseuille_30.xml
  # Poiseuille_60.xml
  # Poiseuille_90.xml
  # Poiseuille_120.xml
  # Poiseuille_180.xml
  # Poiseuille_240.xml

  # TC_10_1.xml
  # TC_20_1.xml
  # TC_30_1.xml
  # TC_40_1.xml
  # TC_50_1.xml

  # TC_10_1.xml
  # TC_20_1.xml
  # TC_30_1.xml
  # TC_40_1.xml
  # TC_50_1.xml
  # CylinderLattice_50.xml
  # CylinderLattice_100.xml
  # CylinderLattice_old_50.xml
  # CylinderLattice_old_100.xml
  
  XML/Cavity/Cavity100_50x50.xml 
  XML/Cavity/Cavity1000_50x50.xml 
  XML/Cavity/Cavity10000_50x50.xml 
  XML/Cavity/Cavity100_100x100.xml 
  XML/Cavity/Cavity1000_100x100.xml 
  XML/Cavity/Cavity10000_100x100.xml 
  XML/Cavity/Cavity100_200x200.xml 
  XML/Cavity/Cavity1000_200x200.xml 
  XML/Cavity/Cavity10000_200x200.xml

)

# Number of processes to use
num_procs=1

# Path to the simulation executable
executable="./sph_dlb"

# Loop over the input files and run the simulation sequentially
for input_file in "${input_files[@]}"; do
    echo "Running simulation with input file: $input_file"
    $executable $input_file

    echo "Simulation for $input_file completed."
done

echo "All simulations completed."