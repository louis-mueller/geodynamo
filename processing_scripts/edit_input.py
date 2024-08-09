# edit_input.py
# -*- coding: utf-8 -*-

'''
Automated Edit of the Input file of CHIC    
Created by Louis Müller (18.04.2024)

Version (09.08.2024)

Summary:
Edits predefined values in list of child directories.
'''

import re

# Input
# --------------------------------------------------------------------
# Define the conditions and their respective values
input_file = 'input.txt'  # Update this path as necessary
parameter = 'M_E'  # You can adjust this based on the condition

conditions = {
    'Condition1': 'NewValue1',
    'Condition2': 'NewValue2',
    'Condition3': 'NewValue3',
}

# Example section and parameter - adjust these as needed
section_name = 'Geometry'
#---------------------------------------------------------------------

def update_parameter(input_file, section_name, parameter, new_value):
    with open(input_file, 'r') as file:
        lines = file.readlines()    # read all lines of input file
    
    updated_lines = []
    section_pattern = re.compile(r'^!\s*' + re.escape(section_name) + r'\s*0.*')
    param_pattern = re.compile(r'^\s*[^!]*\b' + re.escape(parameter) + r'\b\s*=\s*[^!]*')

    in_section = False

    for line in lines:
        if section_pattern.search(line):
            in_section = True
        
        if in_section and param_pattern.search(line):
            line = re.sub(r'^\s*[^!]*\b' + re.escape(parameter) + r'\b\s*=\s*[^!]*', f'{parameter} = {new_value}', line)
            in_section = False
        
        updated_lines.append(line)
    
    with open(input_file, 'w') as file:
        file.writelines(updated_lines)

def main():
    
    for condition, new_value in conditions.items():
        update_parameter(input_file, section_name, parameter, new_value)
    
    print("Input file updated successfully.")

if __name__ == "__main__":
    main()
