import numpy as np
import os
import subprocess
import shutil

base_dir = '/scratch/Simulations_Louis/Debug/'
process_dir = '~/git/geodynamo/processing_scripts/'
process_sh = 'process_debug.sh'
input_file = 'input.txt'
output_file = 'data_prof_M_ 1.0_Fe_ 37.700_FeM_ 10.000_H2O_  0.000.res'
save_dir = '/scratch/Simulations_Louis/Prof/'

Si = np.arange(16.)
M0 = np.arange(21.)

params, name = [], []
for i in Si:
    for j in M0:
        params.append([i, j])
        name.append(f"M1_Fe37-7_wSi0{int(i)}0_M0{int(j)}0.res")

for i in range(len(params)):
    # Read the file and modify lines with specific keywords, preserving alignment
    with open(base_dir + input_file, 'r') as file:
        lines = file.readlines()

    # Modify the file with updated values and write back
    with open(base_dir + input_file, 'w') as file:
        for line in lines:
            stripped_line = line.lstrip()  # Remove leading whitespace for checking keywords
            if stripped_line.startswith('XiS0'):
                # Find the part with 'X_M0', split at '=' and '!' to isolate the value to replace
                parts = line.split('XiS0')
                before = parts[0]
                after = parts[1]
                value_and_comment = after.split('!')
                value_part = value_and_comment[0].split('=')
                
                new_value = f"{params[i][0]:<15}"  # Left-align within 15 characters
                value_and_comment[0] = f"            =       {new_value} "
                
                line = before + 'XiS0' + '!'.join(value_and_comment)

            elif stripped_line.startswith('X_M0'):
                parts = line.split('X_M0')
                before = parts[0]
                after = parts[1]
                value_and_comment = after.split('!')  # Split to isolate value section from comments
                value_part = value_and_comment[0].split('=')  # Split at '=' to locate the numeric value
                
                new_value = f"{params[i][1]:<15}"  # Left-align within 15 characters
                value_and_comment[0] = f"            =       {new_value} "  # Maintain spacing before value
                
                line = before + 'X_M0' + '!'.join(value_and_comment)
            
            # Write the modified line back to the file
            file.write(line)
    print(f"Finished writing parameter set {params[i]} to file.")

    # Construct the full path to the program
    program_path = os.path.join(process_dir, process_sh)

    # Check if the program exists in the directory
    if os.path.isfile(program_path) and os.access(program_path, os.X_OK):
        print(f"{process_sh} found. Running the program...")
        
        # Run the program
        result = subprocess.run([program_path])

        # Check if the program ran successfully
        if result.returncode == 0:
            print("Program ran successfully.")
        else:
            print("Program encountered an error.")
        
        # Check if the output file exists in the base_dir
        output_path = os.path.join(base_dir, output_file)
        if os.path.isfile(output_path):
            print(f"Output file '{output_file}' successfully created in {base_dir}.")
        else:
            print(f"Output file '{output_file}' was not found in {base_dir}.")
        # Move the output file to a new directory
        destination_path = os.path.join(save_dir, name[i])
    
        # Copy and rename the file
        try:
            shutil.copy(output_path, destination_path)
            print(f"Output file '{output_file}' has been successfully copied and renamed to '{name[i]}' in {save_dir}.")
        except Exception as e:
            print(f"An error occurred while copying and renaming the file: {e}")

    else:
        print(f"{process_sh} not found or is not executable in {process_dir}.")