#!/bin/bash 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
# this can be used individually, but it is more effectively being used by the other script file which calls it. 
SUBMITFILE=$1
SRRFILE=$2
samSUBMITFILE=$SUBMITFILE".sam" 
reSUBMITFILE=$SUBMITFILE"-result"
COUNTER=0
# uncomment the reading portions to test out the files next time, since ./reading takes forever 5-10 mins.
g++ reading3.cpp -o reading
wait $!  
# summarySUBMITFILE="summary-"+$SUBMITFILE
./reading $samSUBMITFILE $reSUBMITFILE # running all the files into reading.cpp 
wait $! 
python search_transposeable.py $SUBMITFILE 
wait $! 
# remember the instructions to run phrap 
# shit, if you have your harddrive from school then it's doable here. 
# checking to see if ZI213-breakpoint is empty or not, if not then proceed
wait $! 
targetSUBMITFILE=$SUBMITFILE"-breakpoints"
if [-s $targetSUBMITFILE ]
then 
    GARBAGEFILE="garbagefile"
    # meaning that it only has transposeable elements here.

    CURRENT=$targetSUBMITFILE" is empty and has nothing in it" 
    echo $CURRENT
    $CURRENT>>GARBAGEFILE
else 
    idSUBMITFILE="id"$targetSUBMITFILE".txt"
    # cut f1 > has to occur somewhere here man., 
    cut -f2 < $targetSUBMITFILE > $idSUBMITFILE
    # do the perl and phrap shit here
    fastqSUBMITFILE_1=$SRRFILE"_1.fastq.gz"
    fastqSUBMITFILE_2=$SRRFILE"_2.fastq.gz"
    finishFILE=$SUBMITFILE"_reads.txt"
    gzip -dc $fastqSUBMITFILE_1 | perl parse_reads.pl $idSUBMITFILE > $finishFILE  
    wait $! 
    gzip -dc $fastqSUBMITFILE_2 | perl parse_reads.pl $idSUBMITFILE >> $finishFILE 
    wait $! 
    first_finishFILE=$SUBMITFILE"_1.out"
    second_finishFILE=$SUBMITFILE"_1.out.qual"
    perl fastq2fa_qual.pl < $finishFILE > $first_finishFILE 2> $second_finishFILE
    # start the next part 
    wait $! 
    phrap -vector_bound 0 -forcelevel 10 $first_finishFILE # double check on this file right here. 
    wait $! 
    echo $SUMBITFILE" is done:"
fi 
