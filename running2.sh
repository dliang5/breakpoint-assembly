#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this is called by python because python is alot easier to call
SUBMITFILE=$1
reSUBMITFILE=$SUBMITFILE+"-result"
# summarySUBMITFILE="summary-"+$SUBMITFILE
./reading $SUBMITFILE $reSUBMITFILE # running all the files into reading.cpp 
wait $! 
python search_transposeable.py $SUBMITFILE 
wait $! 
# remember the instructions to run phrap 
# shit, if you have your harddrive from school then it's doable here. 
# checking to see if ZI213-breakpoint is empty or not, if not then proceed
targetSUBMITFILE=$SUBMITFILE+"-breakpoints"
if [-s $targetSUBMITFILE ]
then 
    echo $targetSUBMITFILE + " is empty and has nothing in it" 
else 
    # do the perl and phrap shit here
    fastqSUBMITFILE_1=$SUBMITFILE+"_1.fastq.gz"
    fastqSUBMITFILE_2=$SUBMITFILE+"_2.fastq.gz"
    finishFILE=$SUBMITFILE+"ID.txt" 
    perl parse_reads.pl $fastqSUBMITFILE_1 > $finishFILE 
    perl parse_reads.pl $fastqSUBMITFILE_2 > $finishFILE 
    wait $! 
    # start the next part 
    trashSUBMITFILE="reading"+$SUBMITFILE
    qualSUBMITFILE="reading"+$SUBMITFILE # make sure on this end before submitting this  
    perl fastq2fa_qual.pl $finishFILE > $trashSUBMITFILE 2> $qualSUBMITFILE
    wait $! 
    phrap -vector_bound 10 -force_level 0 <the file name to assemble> # double check on this file right here. 
    wait $! 
    echo $SUMBITFILE + " is done:"
