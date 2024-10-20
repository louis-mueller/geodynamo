## CHIC Processing Scripts
This repository contains processing scripts for CHIC (Noack et al. 2015) input and output files.

Author: Louis Müller (12.08.2024)

### Dependencies for visu_profs.py and process_profs.sh
Unix System Requirements:
- python3
- awk

### Summary vis_main.py
Python visualization suite for core, heat-data, boundary-params and interior structure Profiles
main files to be edited are: vis_main.py and input_data.py
The output file is stored to the ouput directoy (typically: local_dir/Plots/output_file.svg).

Modules: numpy, scipy, matplotlib, paramiko

Necessary module list for curta.zedat.fu-berlin.de to plot locally:
```
 module load Python/3.12.3-GCCcore-13.3.0
```
before plotting you need to set up a virtual enviroment an install the required modules
On a side note CHIC is compiled successfully with intel/2020a on curta so it is usefull to purge your modules 
before executing CHIC otherwise Compiler issues may occur.

### Summary process_profs.sh
Bash Script to automatically run a sequence of interior structure simulations with CHIC. 
The code expects you to define a base directory from where you will run all operations (typically scratch/Simulations_xxx)
The Profiles are then run from the Profile directory which you can define (typically $base_dir/Prof)
In the base directory the program will expect you to have the Simulation directories defined as in conditions. 
If this is not the case it will make the necessary directory named in conditions. 
Make sure you have a working version of the CHIC executable and, 
if you want to use Perple_X, a directory called: "Composition_Final_MORE_FEO" both in your base directory.

### Summary visu_profs.py
Run This program in the base directory as `python3 visu_profs.py "DIR_STRS" "OUTPUT_FILE_NAME" "MAX_COUNT"`
Any MAX_COUNT of directories named as listed in DIR_STRS containing a file named INT_STRUCT_DATA (e.g., profs.res), 
with the data stored in columns, are plotted in a pre defined configuration of Subplots. 
The color is chosen by a letter snippet defining the planet mass (M), and three further snippets can 
be defined in ID to differentiate simulations accordingly.

### Summary process_debug.sh
This script is used to streamline the CHIC testing and debugging workflow.
It can remove previous files in the Debug directory, can compile CHIC, 
copies the executable to the simulation directory, and runs CHIC.

### Setup
The set up code describes process_profs.sh but can be used for other bash scripts as well.
The code was initially developed on a Windows system, which may lead to issues with line endings. To convert line endings from Windows to Unix format, run:
```
dos2unix ./process_profs.sh
```

After converting the line endings, ensure the script has execute permissions with:
```
chmod +x ./process_profs.sh
```

Ensure you have a working version of the CHIC executable. If you are using Perple_X (Connolly 2005), you also need a directory named "Composition_Final_MORE_FEO" (Balduin 2024) containing the phase data in your base directory.
Additionally, download the PREM500.csv file (Durek & Ekström 1996) to the base directory. If PREM is enabled in visu_profs.py, set the variable n to the length of the data.

Edit the directory paths in the Input block of process_profs.sh to match your individual path configuration. 
Customize other input values as needed.

After setting up all directories and configurations, you can execute the bash script with:
```
./process_profs.sh
```

### ToDo
- Add another Input_var (On/OFF) that lets you just plot, and does not create new profs.res files
- Add a color map Input_var that discretizes a contious pallete (flexible coloring).
- Add stable legend positioning when PREM is set to TRUE (flexible Legend).
- Add additional arguments to be passed from the bash script to visu_profs.py.
- Add a better routine for saving previous input.txt files
- Continue testing edge cases for margin use.
