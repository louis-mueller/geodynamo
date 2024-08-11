#!/bin/bash

# Version (11.08.2024), Louis Müller

# Bash Script to automatically run a sequence of interior structure profile simulations with CHIC 
# The code expects you to define a base directory from where you will run all operations (typically scratch/Simulations_xxx)
# The Profiles are then run from the Profile directory which you can define (typically $base_dir/Prof)
# In the base directory the program will expect you to have the Simulation directories defined as in conditions. 
# If this is not the case it will make the necessary directory named in conditions.

# Exit the script on any error
set -e

# Input
#--------------------------------------------------------------------------------------
base_dir="/scratch/Simulations_Louis"
prof_dir="$base_dir/Prof"
chic_executable="$base_dir/CHIC_070824"
output_pattern="$prof_dir/data_prof_M*"
save_old_input=0                                # if greater 0 the edited input is saved as a copy with the correct dir_name
max_count=11                                     # change this value for the amount of simulations you would like to process
count=0
last_size=0                                     # subsequent four definitions are used to check if the output file 
stable_size=0                                   # is still being written into
stable_time=0
max_stable_time=180                             # Maximum time to wait for file stability in seconds

declare -A conditions
conditions=(                                    # conditions (right) of each directory name (left)
    ["M1_Fe30_sFe6-5_p"]="M_E=1.0,X_Fe=30,Dl=100000.0,Dcr0=50000.0"      # Define values for M_E and X_Fe for conditions                                           
    ["M1_Fe60_sFe6-5_p"]="M_E=1.0,X_Fe=60,Dl=100000.0,Dcr0=50000.0"   # The strings are seperated by "_" or "=" and the specific number is read.
    ["M2_Fe30_sFe6-5_p"]="M_E=2.0,X_Fe=30,Dl=72146.0,Dcr0=36073.0"      # in the code: value =$(echo "$condition" | awk -F'[,=]' '{print $2}')
    ["M2_Fe60_sFe6-5_p"]="M_E=2.0,X_Fe=60,Dl=72146.0,Dcr0=36073.0"     # fields (-F) are:     1 = 2 , 3  = 4
    ["M3_Fe30_sFe6-5_p"]="M_E=3.0,X_Fe=30,Dl=59604.0,Dcr0=29802.0"      #                    "M_E=3.0,X_Fe=30"
    ["M3_Fe60_sFe6-5_p"]="M_E=3.0,X_Fe=60,Dl=59604.0,Dcr0=29802.0"     
    ["M4_Fe30_sFe6-5_p"]="M_E=4.0,X_Fe=30,Dl=52051.0,Dcr0=26026.0"
    ["M4_Fe60_sFe6-5_p"]="M_E=4.0,X_Fe=60,Dl=52051.0,Dcr0=26026.0"
    ["M5_Fe30_sFe6-5_p"]="M_E=5.0,X_Fe=30,Dl=46858.0,Dcr0=23429.0"
    ["M5_Fe60_sFe6-5_p"]="M_E=5.0,X_Fe=60,Dl=46858.0,Dcr0=23429.0"
)
#--------------------------------------------------------------------------------------

# Function to check if a variable is an integer or floating-point number and exits if not
check_number() {
    local value="$1" 
    # Regular expression for a valid integer or floating-point number
    if ! [[ "$value" =~ ^-?[0-9]*(\.[0-9]+)?$ ]]; then
        echo "Error: '$value' is not a valid number."
        exit 1
    else 
        echo "$value is a number."
    fi
}

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
        # Loop through each condition
        for dir_name in "${!conditions[@]}"; do
            output_dir="$base_dir/$dir_name"
            echo "$outputdir"

            if [ $count -ge $max_count ]; then 
                echo "Your test count has reached its limit."
                exit 1
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
            echo "Updating input.txt '$dir_name': $condition"

            # Use awk to update input.txt based on the condition
            awk -v condition="${conditions[$dir_name]}" '
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
            
            # --------code segment checks if output file is still being written into-----------------------------
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
            #-----------------------------------------------------------------------------------------------------

            # Check if output file exists and copy it
            if [ -f "$output_file" ]; then
                cp "$output_file" "$output_dir/profs.res" #the Interior Structure data is copied to the correct directory
                echo "Output file copied to $output_dir/profs.res"
            else
                echo "Error: Output file $output_file not found in $prof_dir."
                exit 1
            fi

            echo "Processing for directory $output_dir complete." 

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

echo "All directories processed successfully."