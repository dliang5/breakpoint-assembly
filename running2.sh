#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this is called by python because python is alot easier to call
SUBMITFILE=$1
samSUBMITFILE=$SUBMITFILE+".sam" 
reSUBMITFILE=$SUBMITFILE+"-result"
COUNTER=0
g++ reading3.cpp -o reading 
# summarySUBMITFILE="summary-"+$SUBMITFILE
./reading $samSUBMITFILE $reSUBMITFILE # running all the files into reading.cpp 
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
    idSUBMITFILE="id"+$targetSUBMITFILE+".txt"
    # cut f1 > has to occur somewhere here man., 
    cut -f2 < $targetSUBMITFILE > $idSUBMITFILE
    # do the perl and phrap shit here
    fastqSUBMITFILE_1=$SUBMITFILE+"_1.fastq.gz"
    fastqSUBMITFILE_2=$SUBMITFILE+"_2.fastq.gz"
    finishFILE=$SUBMITFILE+"reads.txt"
    gzip -dc $SRR_1 | perl parse_reads.pl $idSUBMITFILE > $finishFILE  
    wait $! 
    gzip -dc $SRR_2 | perl parse_reads.pl $idSUBMITFILE >> $finishFILE 
    wait $! 
    first_finishFILE=$SUBMITFILE+"_1.out"
    second_finishFILE=$SUBMITFILE+"_1.out.qual"
    perl fastq2fa_qual.pl < $finishFILE > $first_finishFILE 2> $second_finishFILE
    # start the next part 
    wait $! 
    phrap -vector_bound 0 -force_level 10 $first_finishFILE # double check on this file right here. 
    wait $! 
    echo $SUMBITFILE + " is done:"
