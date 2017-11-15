#!/bin/bash
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 

""" This is the offical python script file to run the script file on the entire spreadsheet 
    It reads RAL-README and match the RAL names with the SRR/ERR names in TableS1_individuals.csv with grep

    run : python2.7 officialScript.py 
    output : qsub running3.sh SRR_name 
""" 

import os, subprocess as sp, sys

# reading from john pool's spreadsheet 
# Ex. Name, name1, originalRegion, SRRname, chromosomes

#however, we are going to check RAL files first
with open ("RAL-README", 'w') as f: 
    f.write("this is here to create the RAL file and store all of the corresponding information in there")
    
inputFile = "RAL-README" 
writeFile = open(inputFile, 'w') 

with open("TableS1_individuals.csv", 'r') as poolFile: 
    for line in poolFile: 
        content = [x.strip() for x in line.split(",")]
        setName = content[0] 
        if setName[:3] != 'RAL': 
            # no need for other yet
            continue 
        else: 
            counter = 0 
            checkTwice = 'SRX'
            
            defaultSRRname = 'qsub running3.sh {0} '.format(setName) 
            for position, entry in enumerate(content): 
                if entry[:3] == "SRX": 
                    defaultSRRname = defaultSRRname + "{0} ".format(entry.strip('"'))
            # print (defaultSRRname)
            sp.Popen(defaultSRRname, shell=True, exectuable= '/bin/bash')
            writeFile.write(defaultSRRname + "\n")

print ("done") 
writeFile.close()