#!/bin/bash
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 

""" This is the offical python script file to run the script file on the entire spreadsheet """ 

import os, subprocess as sp, sys

# reading from john pool's spreadsheet 
# Ex. Name, name1, originalRegion, SRRname, chromosomes

#however, we are going to check RAL files first

with open("TableS1_individuals.csv", 'r') as poolFile: 
    sp.Popen("mkdir RAL", shell=True, executable='/bin/bash')
    for line in poolFile: 
        content = [x.strip() for x in line.split(",")]
        if content[:3] != 'RAL': 
            continue 
        else: 
            SRRname = content[3] 
            path = "qsub running3.sh" + SRRname + " >> RAL" 
            sp.Popen(path, shell=True, executable='/bin/bash')
            