#!/bin/bash

# Summary
# This script is used to streamline the CHIC testing and debugging workflow.

# It can remove previous files in the Debug directory, can compile CHIC, 
# copies the executable to the simulation directory, and runs CHIC.

# Exit the script on any error
set -e

# Input
#--------------------------------------------------------------------------------------
base_dir="/scratch/Simulations_Louis"
chic_executable="/home/louismueller/git/chic_v2/bin/CHIC"     # Replace with the actual path to CHIC executable
debug_dir="$base_dir/Debug"                                   # Directory where the simulations will be run should contain input.txt and prof.res (for M1_Fe30_sFe6-5_p)
delete_data_files=1                                           # If greater than 0, files starting with 'data' and a few others will be deleted
make_chic=1                                                   # If greater than 0, CHIC will be compiled

#--------------------------------------------------------------------------------------

echo "Starting CHIC testing and debugging workflow: process_debug.sh..."

# Add make (CHIC) command here
if [[ "$make_chic" -gt 0 ]]; then
    echo "Compiling CHIC..."
    cd /home/louismueller/git/chic_v2 || exit 1
    make
    echo "CHIC compiled successfully."
fi

# Define current date in DDMMYYYY format
current_date=$(date +"%d%m%Y")

# Define the target name for the copied chic executable
target_chic="${base_dir}/CHIC_${current_date}" 

if [ -f "$chic_executable" ]; then
    cp "$chic_executable" "$target_chic"
    echo "CHIC executable copied to ${target_chic}"
else
    echo "Error: CHIC executable not found at specified path."
    exit 1
fi

# Check if Debug directory exists
if [ -d "$debug_dir" ]; then
    echo "Debug directory found. Checking for files starting with 'data'..."

    data_files=$(find "$debug_dir" -maxdepth 1 -type f -name 'data*')
    fort_files=$(find "$debug_dir" -maxdepth 1 -type f -name 'fort*')
    result_files=$(find "$debug_dir" -maxdepth 1 -type f -name 'all_results**')
    
    # Remove all previous data files
    if [[ -n "$data_files" && "$delete_data_files" -gt 0 ]]; then

        cd "$debug_dir" || exit 1

        if [[ -n "$fort_files" ]]; then
            echo "Removing other irrelevant non-datafiles:"
            rm -r ./fort* 
        fi

        if [[ -n "$result_files" ]]; then
            echo "Removing other irrelevant non-datafiles:"
            rm -r ./all_results* 
        fi

        echo "Removing files starting with 'data' in $(pwd):"
    
        while IFS= read -r file; do
            echo "$file"
            rm "$file"
        done <<< "$data_files"
        echo "Files removed successfully."
    else
        echo "No files starting with 'data' found in $debug_dir."
    fi
else
    echo "Error: The dirctory: $debug_dir does not exist."
    exit 1
fi

# Continue by running chic code in the debug directory
echo "Running CHIC in the Debug directory..."

# Change to Debug directory
cd "$debug_dir" || exit 1
echo "The current working directory is: $(pwd)"

# Run CHIC program
"$target_chic"

echo "CHIC run completed."