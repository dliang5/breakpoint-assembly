#!/bin/bash
#SBATCH --job-name -mosquito-job
#SBATCH --nodes 1 
#SBATCH --output mosquitotest.out
#SBATCH --ntasks=3

# this calls the other five files to read and parse the files for assembly 
# after this then we are done for now 