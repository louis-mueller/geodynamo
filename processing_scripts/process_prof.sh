#!/bin/bash

# Processing of Interior Structure Simulations 
# Created by Louis Müller (11.08.2024)

# Version (12.08.2024)

# Summary:
# Bash Script to automatically run a sequence of interior structure simulations with CHIC. 
# The code expects you to define a base directory from where you will run all operations (typically scratch/Simulations_xxx)
# The Profiles are then run from the Profile directory which you can define (typically $base_dir/Prof)
# In the base directory the program will expect you to have the Simulation directories defined as in conditions. 
# If this is not the case it will make the necessary directory named in conditions.

# Dependancies: python3, os, sys, matplotlib, numpy, and scipy.
# Make sure you have a working version of the CHIC executable and, 
# if you want to use Perple_X, a directory called: 
# "Composition_Final_MORE_FEO" both in your base directory.

# Exit the script on any error
set -e

# Input
#--------------------------------------------------------------------------------------
base_dir="/scratch/Simulations_Louis"
prof_dir="$base_dir/Prof"
chic_executable="$base_dir/CHIC_070824"
output_pattern="$prof_dir/data_prof_M*"
save_old_input=0                                    # if greater 0 the edited input is saved as a copy with the correct dir_name
max_count=10                                        # change this value for the amount of simulations you would like to process
plot_output=1                                       # if greater 0 plotting (eg. visu_profs.py) code is run as well
plot_dir="/home/louismueller/bin/visu_profs.py"     # path to plotting code 
OUTPUT_FILE_NAME="all_sim_IntStruct_120824.svg"     # passed to plotting code 

dir_strs=(
    "M1_Fe30_sFe6-5_p" "M1_Fe60_sFe6-5_p" "M2_Fe30_sFe6-5_p" "M2_Fe60_sFe6-5_p"
    "M3_Fe30_sFe6-5_p" "M3_Fe60_sFe6-5_p" "M4_Fe30_sFe6-5_p" "M4_Fe60_sFe6-5_p"
    "M5_Fe30_sFe6-5_p" "M5_Fe60_sFe6-5_p"
)                                                   # dir_strs are passed to plotting code later as well 

param_cond=(
    "M_E=1.0,X_Fe=30,Dl=100000.0,Dcr0=50000.0" "M_E=1.0,X_Fe=60,Dl=100000.0,Dcr0=50000.0"
    "M_E=2.0,X_Fe=30,Dl=72146.0,Dcr0=36073.0" "M_E=2.0,X_Fe=60,Dl=72146.0,Dcr0=36073.0"
    "M_E=3.0,X_Fe=30,Dl=59604.0,Dcr0=29802.0" "M_E=3.0,X_Fe=60,Dl=59604.0,Dcr0=29802.0"
    "M_E=4.0,X_Fe=30,Dl=52051.0,Dcr0=26026.0" "M_E=4.0,X_Fe=60,Dl=52051.0,Dcr0=26026.0"
    "M_E=5.0,X_Fe=30,Dl=46858.0,Dcr0=23429.0" "M_E=5.0,X_Fe=60,Dl=46858.0,Dcr0=23429.0" 
)
#----------------------------------------------------------------------------------------

# Check if the number of directories matches the number of conditions
if [ ${#dir_strs[@]} -ne ${#param_cond[@]} ]; then
    echo "Error: The number of directory strings does not match the number of parameter conditions."
    exit 1
fi

# Create an associative array to hold the key-value pairs
declare -A conditions

# Populate the associative array
for i in "${!dir_strs[@]}"; do
    dir="${dir_strs[$i]}"
    param="${param_cond[$i]}"
    conditions["$dir"]="$param"
done

# Check if Prof directory exists
if [ -d "$prof_dir" ]; then
    echo "Prof directory found. Checking for files starting with 'data'..."

    data_files=$(find "$prof_dir" -maxdepth 1 -type f -name 'data*')

    if [ -n "$data_files" ]; then
        echo "Removing files starting with 'data':"
        # Process each file individually
        while IFS= read -r file; do
            echo "$file"
            rm "$file"
        done <<< "$data_files"
        echo "Files removed successfully."
    else
        echo "No files starting with 'data' found in $prof_dir."
    fi
else
    echo "Error: The dirctory ($prof_dir) does not exist."
    exit 1
fi

# Check if input file exists in Prof directory
if [ -f "$prof_dir/input.txt" ]; then
    echo "Input file found. Processing..."

    # Check if CHIC executable exists in the base directory
    if [ -x "$chic_executable" ]; then

        DIR_STRS=""
        count=0
        # Loop through each condition
        for dir_name in "${dir_strs[@]}"; do
            
            data_files=$(find "$prof_dir" -maxdepth 1 -type f -name 'data*')

            if [ -n "$data_files" ]; then
                # Process each file individually
                while IFS= read -r file; do
                    rm "$file"
                done <<< "$data_files"
            fi

            param="${conditions[$dir_name]}"
            output_dir="$base_dir/$dir_name"
            echo "$outputdir"

            if [ $count -ge $max_count ]; then 
                echo "Your test count has reached its limit."
                break
            fi

            # Check if output directory exists
            if [ ! -d "$output_dir" ]; then
                echo "Creating output directory: $output_dir"
                mkdir -p "$output_dir"
            else
                echo "$output_dir exists in the base directory."
            fi
            
            echo "Processing directory: $output_dir"

            # Save the old input.txt file version (optional)
            if [ $save_old_input -gt 0 ]; then
                cp "$prof_dir/input.txt" "$prof_dir/input_IntStruct_before.txt"
                echo "Earlier input.txt file saved as input_IntStruct_before.txt"
            fi
            
            # Modify input.txt file based on the current condition
            echo "Updating input.txt '$dir_name': '$param'"

            # Use awk to update input.txt based on the condition
            awk -v condition=$param '
            BEGIN {
                # Initialize field separator and output field separator
                FS="[ \t]*=[ \t]*";   # Handle possible spaces and tabs around "="
                OFS="\t\t\t=\t\t";    # Set specific output structure 

                # Split the condition into parameter-value pairs
                split(condition, pairs, ",");
                for (i in pairs) {
                    split(pairs[i], kv, "=");
                    params[kv[1]] = kv[2];
                }
            }
            {
                # Check and update fields based on the parameters
                if ($1 in params) {
                    $2 = params[$1];
                }
                print $0;
            }
            ' "$prof_dir/input.txt" > "$prof_dir/input_temp.txt" && mv "$prof_dir/input_temp.txt" "$prof_dir/input.txt"

            echo "finished editing input.txt for $dir_name"
            
            # Change to Prof directory
            cd "$prof_dir" || exit

            # Run CHIC program
            "$chic_executable"

            # Find the output file matching the pattern (assuming there is only one)
            output_file=$(find "$prof_dir" -maxdepth 1 -type f -name 'data_prof_M*' -print -quit)

            # Check if the output file was found
            if [ -z "$output_file" ]; then
                echo "No matching files found."
                exit 1
            fi
            echo "Found output file: $output_file"

            # Check if output file exists and copy it
            if [ -f "$output_file" ]; then
                cp "$output_file" "$output_dir/profs.res" #the Interior Structure data is copied to the correct directory
                echo "Output file copied to $output_dir/profs.res"
            else
                echo "Error: Output file $output_file not found in $prof_dir."
                exit 1
            fi

            echo "Processing for directory $output_dir complete." 

            # Check if plotting is enabled
            if [ $plot_output -gt 0 ]; then
                # If DIR_STRS is not empty, append a comma before adding the new value
                if [ -n "$DIR_STRS" ]; then
                    DIR_STRS+=","
                fi
                # Append the current dir_str to DIR_STRS
                DIR_STRS+="$dir_name"
            fi

            count=$(($count + 1)) # counts the amount of directories processed
        
        done

        echo "CHIC executed successfully $count time*s."

        # Navigate back to the base directory
        cd "$base_dir" || exit
    else
        echo "Error: CHIC executable not found in $base_dir."
        exit 1
    fi
else
    echo "Error: input.txt not found in $prof_dir."
    exit 1
fi

if [ $plot_output -gt 0 ]; then
    echo "Plotting started..."

    # Check if input file exists in Prof directory
    if [ -f "$plot_dir" ]; then
        echo "Python plotting script found. Processing..."

        MAX_COUNT=$max_count

        echo "Directories passed to plotting script: $DIR_STRS"
        echo "Output file name passed to plotting script: $OUTPUT_FILE_NAME"
        echo "Upper limit of processed Directories passed to plotting script: $MAX_COUNT"

        # Call the Python script and passes the correct arguments
        python3 $plot_dir "$DIR_STRS" "$OUTPUT_FILE_NAME" "$MAX_COUNT"
    else 
        echo "Error: File at $plot_dir was not found"
        exit 1
    fi

fi

echo "All directories processed successfully." 