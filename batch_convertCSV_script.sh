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
    # Poiseuille_old_summation_12_30_1rf_1prc
    # Poiseuille_old_summation_24_60_1rf_1prc
    # Poiseuille_old_summation_36_90_1rf_1prc
    # Poiseuille_old_summation_48_120_1rf_1prc
    # Poiseuille_old_summation_72_180_1rf_1prc
    # Poiseuille_old_summation_96_240_1rf_1prc

    # Poiseuille_new_summation_12_30_1rf_1prc
    # Poiseuille_new_summation_24_60_1rf_1prc
    # Poiseuille_new_summation_36_90_1rf_1prc
    # Poiseuille_new_summation_48_120_1rf_1prc
    # Poiseuille_new_summation_72_180_1rf_1prc
    # Poiseuille_new_summation_96_240_1rf_1prc
    
    Poiseuille_old_differential_12_30_1rf_1prc
    Poiseuille_old_differential_24_60_1rf_1prc
    Poiseuille_old_differential_36_90_1rf_1prc
    Poiseuille_old_differential_48_120_1rf_1prc
    Poiseuille_old_differential_72_180_1rf_1prc
    Poiseuille_old_differential_96_240_1rf_1prc

    Poiseuille_new_differential_12_30_1rf_1prc
    Poiseuille_new_differential_24_60_1rf_1prc
    Poiseuille_new_differential_36_90_1rf_1prc
    Poiseuille_new_differential_48_120_1rf_1prc
    Poiseuille_new_differential_72_180_1rf_1prc
    Poiseuille_new_differential_96_240_1rf_1prc

    # Poiseuille_new_summation_12_30_2rf_1prc
    # Poiseuille_new_summation_24_60_2rf_1prc
    # Poiseuille_new_summation_36_90_2rf_1prc
    # Poiseuille_new_summation_48_120_2rf_1prc
    # Poiseuille_new_summation_72_180_2rf_1prc
    # Poiseuille_new_summation_96_240_2rf_1prc

    # TaylorCouette_new_summation_10_10_1rf_1prc
    # TaylorCouette_new_summation_20_20_1rf_1prc
    # TaylorCouette_new_summation_30_30_1rf_1prc
    # TaylorCouette_new_summation_40_40_1rf_1prc

    # TaylorCouette_new_summation_10_10_2rf_1prc
    # TaylorCouette_new_summation_20_20_2rf_1prc
    # TaylorCouette_new_summation_30_30_2rf_1prc
    # TaylorCouette_new_summation_40_40_2rf_1prc

    # TaylorCouette_old_summation_10_10_1rf_1prc
    # TaylorCouette_old_summation_20_20_1rf_1prc
    # TaylorCouette_old_summation_30_30_1rf_1prc
    # TaylorCouette_old_summation_40_40_1rf_1prc


####################################################################33
    # TaylorCouette_new_differential_10_10_1rf_1prc
    # TaylorCouette_new_differential_20_20_1rf_1prc
    # TaylorCouette_new_differential_30_30_1rf_1prc
    # TaylorCouette_new_differential_40_40_1rf_1prc
    # TaylorCouette_new_differential_50_50_1rf_1prc

    # TaylorCouette_new_differential_10_10_2rf_1prc
    # TaylorCouette_new_differential_20_20_2rf_1prc
    # TaylorCouette_new_differential_30_30_2rf_1prc
    # TaylorCouette_new_differential_40_40_2rf_1prc
    # TaylorCouette_new_differential_50_50_2rf_1prc
    
    # TaylorCouette_new_summation_10_10_1rf_1prc
    # TaylorCouette_new_summation_20_20_1rf_1prc
    # TaylorCouette_new_summation_30_30_1rf_1prc
    # TaylorCouette_new_summation_40_40_1rf_1prc
    # TaylorCouette_new_summation_50_50_1rf_1prc


    # TaylorCouette_old_differential_10_10_1rf_1prc
    # TaylorCouette_old_differential_20_20_1rf_1prc
    # TaylorCouette_old_differential_30_30_1rf_1prc
    # TaylorCouette_old_differential_40_40_1rf_1prc
    # TaylorCouette_old_differential_50_50_1rf_1prc

    # TaylorCouette_new_differential_10_10_0rf_1prc
    # TaylorCouette_new_differential_20_20_0rf_1prc
    # TaylorCouette_new_differential_30_30_0rf_1prc
    # TaylorCouette_new_differential_40_40_0rf_1prc
    # TaylorCouette_new_differential_50_50_0rf_1prc

    # TaylorCouette_new_differential_flatvolume_10_10_1rf_1prc
    # TaylorCouette_new_differential_flatvolume_20_20_1rf_1prc
    # TaylorCouette_new_differential_flatvolume_30_30_1rf_1prc
    # TaylorCouette_new_differential_flatvolume_40_40_1rf_1prc
    # TaylorCouette_new_differential_flatvolume_50_50_1rf_1prc


) 

# Loop through the hardcoded list of files
for filename in "${files[@]}"; do
    # Call the Python script for each file pattern
    echo "Processing files for pattern: $filename"
    pvpython convertCSV.py "$filename" $start $end
done