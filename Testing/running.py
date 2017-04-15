#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
"""
    This is the head file that will submit running2.sh individually with 1 arguments
    python running.py <set_number (my bad, ill change it in the future) > 

    basically, it will read each line from the set_*.txt to get the SRRname when using perl to 
    extract the SRR_IDs to be assemble. 
"""

import os, subprocess as sp, sys, glob
print 'argument list:', str(sys.argv) 
arg_content = str(sys.argv).split()
# name of the sam file to be read into running2.sh and the corresponding sites
number = arg_content[1].lstrip("['").rstrip(",']")
print (number) 
name = "set_"+number+"_"
filename = "set_"+number+".txt"
lineNumber=0 
# read the SRRfile and the content inside. 
with open(filename, 'r') as f: 
    # lineNumber = sum(1 for _ in f) 
    # f.seek(0) 
    for count, j in enumerate(f): # getting the SRRname here 
        """ 
        here, this then creates the submission to qsub to the cluster 
        so it'll be "qsub running2.sh <set_1_1.sam> <SRRname of set_1_1.sam"> 
        the SRRname comes from the set_*.txt 
        """
        count+=1
        submission = name+str(count) # this forms set_1_* <* is any number> 
        path = "qsub running2.sh " + submission + " " + j.strip("\n") 
        s_path = '/bin/bash'
        print(submission + " ->" +  j.strip("\n") )  
        sp.Popen(path, shell=True,executable=s_path ) 
        # print ("set_"+"1_"+str(count)+" to -> " + i).strip("\n")

# for i in xrange(0,20):
#     counter = i+1
#     submission = "set_"+number+str(counter)+".sam"
#     path = "qsub running2.sh " + submission # qsub running2.sh set_1_1.sam etc  
#     subprocess.Popen(path)
