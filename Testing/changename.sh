#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this is just to change the name of all the set files into their respective SRR name 

##### -> WARNING : i have not tested this file out so try it out before going at it. 
SET=$2 # what set number this is
rm running.sh.o* 
rm "summary-set_"$SET_"*-*"
rm "set_"$SET"_*_reads.txt"
rm garbage_set_* 
rm good_set_*