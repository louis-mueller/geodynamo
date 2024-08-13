## CHIC Processing Scripts
This repository contains processing scripts for CHIC (Noack et al. 2015) input and output files.

Author: Louis Müller (12.08.2024)

### Dependencies for visu_profs.py and process_profs.sh
Unix System Requirements:
- python3
- awk

Python Environment Requirements:
- os
- sys
- matplotlib
- numpy
- scipy

### Setup
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
- Add stable legend positioning when PREM is set to TRUE.
- Add additional arguments to be passed from the bash script to visu_profs.py.
- Continue testing edge cases for margin use.
- Add a better routine for saving previous input.txt files 
