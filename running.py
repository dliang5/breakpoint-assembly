#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this is to run reading3.cpp, search_transposeable.py, and phrap in one go
import os
import subprocess
import sys 
print 'argument list:', str(sys.argv) 
arg_content = str(sys.argv).split()
number = arg_content[1].lstrip("'").rstrip("']") 
for i in xrange(0,20):
    counter = i+1
    submission = "set_"+number+str(counter)+".sam"
    path = "qsub running2.sh " + submission # qsub running2.sh set_1_1.sam etc  
    subprocess.Popen(path)