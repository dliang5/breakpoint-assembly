#!/bin/bash
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 

""" This is the offical python script file to run the script file on the entire spreadsheet """ 

import os, subprocess as sp, sys

# reading from john pool's spreadsheet 
# Ex. Name, name1, originalRegion, SRRname, chromosomes

#however, we are going to check RAL files first
with open ("RAL-README", 'w') as f: 
    f.write("this is here to create the RAL file and store all of the corresponding information in there")
    
inputFile = "RAL-README" 
writeFile = open(inputFile, 'a') 

with open("TableS1_individuals.csv", 'r') as poolFile: 
    sp.Popen("mkdir RAL", shell=True, executable='/bin/bash')
    for line in poolFile: 
        content = [x.strip() for x in line.split(",")]
        setName = content[0] 
        if setName[:3] != 'RAL': 
            # no need for other yet
            continue 
        else: 
            checkTwice = 'SR'
            srrTwice = content[4]  
            if ( srrTwice[:2] == checkTwice): 
                SRRname = content[3]
                SRRname1 = content[4] 
                print (setName,content[3], content[4])
                path = "qsub running3.sh " + setName + " " + SRRname + " " + SRRname1
                sp.Popen(path, shell=True, executable= '/bin/bash')
                writeFile.write("{0} {1} {2}\n".format(setName, SRRname, SRRname1))
                    
            else: 
                SRRname = content[3] 
                path = "qsub running3.sh " + setName+ " " +SRRname  + " " + 0
                # goes to the first one
                print (setName, content[3])
                sp.Popen(path, shell=True, executable='/bin/bash')
                writeFile.write("{0} {1}\n".format(setName, SRRname)) 

writeFile.close()