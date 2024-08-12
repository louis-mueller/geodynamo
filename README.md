# This is a collection of processing scripts for CHIC (Noack et al. 2015) input and output files.

The code is written and updated so far exclusivly by Louis Müller (12.08.2024)

Dependancies for visu_profs.py and process_profs.sh: python3, os, sys, matplotlib, numpy, and scipy.
Make sure you have a working version of the CHIC executable and, if you want to use Perple_X (Connolly 2005), 
a directory called: "Composition_Final_MORE_FEO" (Balduin 2024) containing the phase data both in your base directory.

ToDo:
- Add a stable Legend positioning when PREM is set to TRUE
- Add further args to be passed form the bash script to visu_profs.py
- Continue limit testing for margin use-cases
