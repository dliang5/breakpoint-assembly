# this is the third version of running2.sh which uses officialScript.py 
#!/bin/bash 

# this should do all of the processes from downloading the SRR file to running the assembly.
srrName=$2
srrName1=$3 # this is optional variable
setName=$1
myArgs=( "$@" ) # this holds all of the command line arguments, look at commandLineArgs.sh for more info
setFile="${setName}.bam"

fullName=""
# echo "${myArgs[3]}"
# getting the full name to eventually merge them together later on. 
for i in "${!myArgs[@]}"; do
    if [$i -ne 0] 
    then 
        fullName="${fullName} $i.bam" 
        # doing the process here of aligning and mapping/sorting
        fastq1="${myArgs[i]}_1.fastq.gz"
        fastq2="${myArgs[i]}_2.fastq.gz" 
        bamName="${myArgs[i]}.bam" 
        sortedBam="${myArgs[i]}_sorted.bam" 
        fastq-dump --gzip --split-files "${myArgs[i]}"
        bwa mem ../../DGRP/1_SRR_file/chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $bamName 
    fi
done

if ["${!myArgs[@]}" -ne 1] then 
    samtools merge "${myArgs[0]}.bam" $fullName
    samtools view -h -o "${myArgs[0]}.sam" "${myArgs[0]}.bam" 
    ./reading "${myArgs[0]}.sam" "${myArgs[0]}-result"
    python search_transposeable.py "${myArgs[0]}"
else
    for i in "${!myArgs[@]}"; do
        if [$i -ne 0]
        then 
            samtools view -h -o "${myArgs[0]}.sam" "${myArgs[0]}.bam"
            ./reading "${myArgs[0]}.sam" "${myArgs[0]}-result"
            python search_transposeable.py "${myArgs[0]}"
fi 