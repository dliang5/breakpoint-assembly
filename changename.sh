#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this is just to change the name of all the set files into their respective SRR name 

##### -> WARNING : i have not tested this file out so try it out before going at it. 
SET=$2 # what set number this is
COUNTER=1 # what the counter of the file is   
while IFS ='' read -r line || [[ -n "$line" ]]; do 
    # reading into the file here m8
    new_name="set_"+$SET+"_"+$COUNTER 
    mv $new_name $line # changing the name here 
    let COUNTER=COUNTER+1 # incrementing the counter here.
done < "$1"
