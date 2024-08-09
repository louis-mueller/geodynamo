#!/bin/bash

#SBATCH --job-name=hot_test_1        
#SBATCH --mail-user=louim95@zedat.fu-berlin.de   
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=250                      # not quite sure probably alot smaller
#SBATCH --time=96:00:00                         # each chic sim should take around a day or two
#SBATCH --qos=standard                          # replace with value for your job

#module load GCC
module load intel/2022a

../CHIC                          


