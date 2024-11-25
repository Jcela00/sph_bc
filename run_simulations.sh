#!/bin/bash

# Array of input files
input_files=(
XML/Triangle/Re100/CustomTriangle25.xml
XML/Triangle/Re1000/CustomTriangle25.xml
XML/Ellipse/Re100/CustomEllipse25.xml
XML/Ellipse/Re1000/CustomEllipse25.xml
XML/Triangle/Re100/CustomTriangle50.xml
XML/Triangle/Re1000/CustomTriangle50.xml

# XML/Ellipse/Re100/CustomEllipse50.xml
# XML/Ellipse/Re1000/CustomEllipse50.xml
# XML/Triangle/Re1000/CustomTriangle100.xml
# XML/Ellipse/Re1000/CustomEllipse100.xml
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