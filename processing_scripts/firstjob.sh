#!/bin/bash

#SBATCH --job-name=M1_Fe30_Fe6-5_p        
#SBATCH --mail-user=louim95@zedat.fu-berlin.de   
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=350                       # not quite sure probably alot smaller
#SBATCH --time=96:00:00                         # each chic sim should take around a day or two
#SBATCH --qos=standard                          # replace with value for your job

#module load GCC
module load intel/2022a

cd /scratch/louim95/

# Finds the file that contains CHIC
CHIC=$(find /scratch/louim95/ -name "CHIC")

if [ -z "$CHIC" ]; then
    echo "CHIC not found"
    exit 1
else
    cd /scratch/louim95/M1_Fe30_Fe6-5_p
    $CHIC
fi