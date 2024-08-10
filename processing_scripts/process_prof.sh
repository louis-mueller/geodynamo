#!/bin/bash

# Exit the script on any error
set -e

# Set your working directory and paths
base_dir="/scratch/Simulations_Louis"
prof_dir="$base_dir/Prof"
chic_executable="$base_dir/CHIC_070824"
output_pattern="$prof_dir/data_prof_M*"
count=0
max_count=10 # change this value for the amount of simulations you would like to process


# Define values for M_E and X_Fe for conditions
# The strings are seperated by _ or = and the specific number is read.
declare -A conditions
conditions=(
    ["M1_Fe30_sFe6-5_p"]="M_E=1_X_Fe=30"
    ["M1_Fe60_sFe6-5_p"]="M_E=1_X_Fe=60"
    ["M2_Fe30_sFe6-5_p"]="M_E=2_X_Fe=30"
    ["M2_Fe60_sFe6-5_p"]="M_E=2_X_Fe=60"
    ["M3_Fe30_sFe6-5_p"]="M_E=3_X_Fe=30"
    ["M3_Fe60_sFe6-5_p"]="M_E=3_X_Fe=60"
    ["M4_Fe30_sFe6-5_p"]="M_E=4_X_Fe=30"
    ["M4_Fe60_sFe6-5_p"]="M_E=4_X_Fe=60"
    ["M5_Fe30_sFe6-5_p"]="M_E=5_X_Fe=30"
    ["M5_Fe60_sFe6-5_p"]="M_E=5_X_Fe=60"
)

# Function to check if a value in string is an integer
is_integer() {
    [[ "$1" =~ ^-?[0-9]+$ ]]
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
    echo "Error: Prof directory does not exist."
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

            if [ $count -gt 0 ]; then 
                echo "our test count has reached its limit."
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
            
            # Modify input.txt file based on the current condition
            condition="${conditions[$dir_name]}"

            m_e_val=$(echo "$condition" | awk -F'[_=]' '{print $3}')
            x_fe_val=$(echo "$condition" | awk -F'[_=]' '{print $6}')

            # Ensure that m_e_val and x_fe_val are integers
            if ! is_integer "$m_e_val"; then
                echo "Error: M_E value '$m_e_val' is not a valid integer."
                exit 1
            fi

            if ! is_integer "$x_fe_val"; then
                echo "Error: X_Fe value '$x_fe_val' is not a valid integer."
                exit 1
            fi
            
            echo "Updating input.txt with M_E=$m_e_val and X_Fe=$x_fe_val..."
            awk -v m_e_val="$m_e_val" -v x_fe_val="$x_fe_val" '
            BEGIN {
                FS = "=";
                OFS = "=";
            }
            # Match lines with M_E and X_Fe, allowing for any amount of whitespace
            /^\s*M_E\s*=/ { $2 = m_e_val }
            /^\s*X_Fe\s*=/ { $2 = x_fe_val }
            { print }
            ' "$prof_dir/input.txt" > "$prof_dir/input_temp.txt" && mv "$prof_dir/input_temp.txt" "$prof_dir/input.txt"

            # Save the modified input.txt file (optional)
            cp "$prof_dir/input.txt" "$prof_dir/input_IntStruct_$dir_name.txt"
            
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
