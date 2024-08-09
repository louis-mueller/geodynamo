#!/bin/bash

# Exit the script on any error
set -e

# Set your working directory and paths
base_dir="$(dirname "$(realpath "$0")")" #this should /scratch/Simulatios_Louis in Thor
prof_dir="$base_dir/Prof"
chic_dir="$base_dir/../CHIC"
output_dir_base="$base_dir/../Output"

# Define preconditions and iterations
preconditions_count=3   # Number of predefined conditions
condition_array=("Condition1" "Condition2" "Condition3")

# Check if Prof directory exists
if [ ! -d "$prof_dir" ]; then
    echo "Error: Prof directory does not exist."
    exit 1
fi

# Process each condition
for condition in "${condition_array[@]}"; do
    echo "Processing condition: $condition"

    # Check if input file exists in Prof directory
    if [ -f "$prof_dir/input.txt" ]; then
        echo "Input file found. Processing..."

        # Run CHIC program
        cd "$chic_dir" || exit
        if [ -f "CHIC" ]; then
            ./CHIC
        else
            echo "Error: CHIC executable not found in $chic_dir."
            exit 1
        fi

        # Wait for output file to be created and fully written
        output_file="$chic_dir/data_prof_M(...).txt"
        last_size=0
        stable_size=0
        stable_time=0
        max_stable_time=180  # Maximum time to wait for file stability in seconds

        while true; do
            current_size=$(stat -c%s "$output_file" 2>/dev/null || echo 0)
            if [ "$current_size" -eq "$last_size" ]; then
                stable_size=$((stable_size + 1))
            else
                stable_size=0
            fi

            if [ "$stable_size" -gt 5 ]; then
                echo "$output_file has stabilized."
                break
            fi

            if [ "$stable_time" -ge "$max_stable_time" ]; then
                echo "Error: $output_file did not stabilize within $max_stable_time seconds."
                exit 1
            fi

            last_size=$current_size
            sleep 10
            stable_time=$((stable_time + 10))
        done

        # Check if output file exists and copy it
        if [ -f "$output_file" ]; then
            output_dir="$output_dir_base/$condition"
            mkdir -p "$output_dir"
            cp "$output_file" "$output_dir/profs.res"
        else
            echo "Error: Output file $output_file not found in $chic_dir."
            exit 1
        fi

        # Delete the original data file in Prof directory
        rm -f "$prof_dir/data_prof_M(...).txt"

        # Modify input.txt file based on precondition
        echo "Modifying input.txt file with condition $condition..."
        # (Here you need to include your specific modification commands for input.txt)
        echo "$condition" > "$prof_dir/input.txt"

        # Navigate back to the base directory
        cd "$base_dir" || exit
    else
        echo "Error: input.txt not found in $prof_dir."
        exit 1
    fi
done

echo "Processing complete."
