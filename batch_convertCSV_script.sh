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
    # Poiseuille_new_differential_ConstV_12_30_1rf_1prc   
    # Poiseuille_new_differential_ConstV_24_60_1rf_1prc
    # Poiseuille_new_differential_ConstV_36_90_1rf_1prc
    # Poiseuille_new_differential_ConstV_48_120_1rf_1prc
    # Poiseuille_new_differential_ConstV_72_180_1rf_1prc
    # Poiseuille_new_differential_ConstV_96_240_1rf_1prc   

    # TaylorCouette_new_differential_constV_10_10_1rf_1prc
    # TaylorCouette_new_differential_constV_20_20_1rf_1prc
    # TaylorCouette_new_differential_constV_30_30_1rf_1prc
    # TaylorCouette_new_differential_constV_40_40_1rf_1prc
    # TaylorCouette_new_differential_constV_50_50_1rf_1prc

    # CylinderLattice_new_summation_50_50_1rf_1prc
    # probes_0_CylinderLattice_new_summation_50_50_1rf_1prc
    # probes_1_CylinderLattice_new_summation_50_50_1rf_1prc

    # CylinderLattice_new_summation_100_100_1rf_1prc
    # probes_0_CylinderLattice_new_summation_100_100_1rf_1prc
    # probes_1_CylinderLattice_new_summation_100_100_1rf_1prc

    # Cavity_new_summation_Re100_50_50_1rf_1prc
    # # probes_0_Cavity_new_summation_Re100_50_50_1rf_1prc
    # # probes_1_Cavity_new_summation_Re100_50_50_1rf_1prc

    # Cavity_new_summation_Re100_100_100_1rf_1prc
    # # probes_0_Cavity_new_summation_Re100_100_100_1rf_1prc
    # # probes_1_Cavity_new_summation_Re100_100_100_1rf_1prc

    # Cavity_new_summation_Re100_200_200_1rf_1prc
    # # probes_0_Cavity_new_summation_Re100_200_200_1rf_1prc
    # # probes_1_Cavity_new_summation_Re100_200_200_1rf_1prc

    # Cavity_new_summation_Re1000_50_50_1rf_1prc
    # # probes_0_Cavity_new_summation_Re1000_50_50_1rf_1prc
    # # probes_1_Cavity_new_summation_Re1000_50_50_1rf_1prc

    # Cavity_new_summation_Re1000_100_100_1rf_1prc
    # probes_0_Cavity_new_summation_Re1000_100_100_1rf_1prc
    # probes_1_Cavity_new_summation_Re1000_100_100_1rf_1prc

    # Cavity_new_summation_Re1000_200_200_1rf_1prc
    # probes_0_Cavity_new_summation_Re1000_200_200_1rf_1prc
    # probes_1_Cavity_new_summation_Re1000_200_200_1rf_1prc

    Cavity_new_summation_Re10000_200_200_1rf_1prc
    # probes_0_Cavity_new_summation_Re10000_200_200_1rf_1prc
    # probes_1_Cavity_new_summation_Re10000_200_200_1rf_1prc

) 

# Loop through the hardcoded list of files
for filename in "${files[@]}"; do
    # Call the Python script for each file pattern
    echo "Processing files for pattern: $filename"
    pvpython convertCSV.py "$filename" $start $end
done