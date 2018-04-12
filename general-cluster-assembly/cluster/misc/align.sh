#!/bin/bash
#SBATCH--job-name=script 
#SBATCH--ntasks=3
#SBATCH--output=file_name.out
# run : sbatch align.sh readfile 
# this is just to download and align files that's all up to bam file.
# it reads in a list of the SRR/ERR names from Drosophila Nexus 

:' DISCLAIMER: this assumes bwa index <ref.fa> has been ran
               <chrom_ref> is the output of bwa index <ref.fa> 
'

# reads the file in the command line args 
while IFS='' read -r line|| [[ -n "$line" ]]; do
    # echo "text from file " $line
    a=( $line )
    # {a[0]} = name of file and {a[1...n]} = SRRnames and erx ones
    srrNames1="${a[0]}"
    srrNames=""
    mergedsrrNames=""
    for word in $line; do
        # echo "$word"
        if [[ $word == SRR* ]] || [[ $word == ERR* ]]
        then
            srrNames+=" ${word}"
            fastq1="${word}_1.fastq.gz"
            fastq2="${word}_2.fastq.gz"
            bamName="${word}.bam"
            sortedBam="${word}_sorted.bam"
            mergedsrrNames+="${bamName} "
            fastq-dump --gzip --split-files $word
            wait $!
            bwa mem chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $bamName
            wait $!
        fi
    done
    b=( $srrNames )
    # this is for multiple SRRnames in one line as seen in drosophila nexus or the RAL-file
    if [ ${#a[@]} -ne 2 ] 
    then
        mergedName=${a[0]}".bam"
        mergedSamName=${a[0]}".sam"
        sortedMergedName=${a[0]}"_sorted.bam"
        resultMergedRal=${a[0]}"-result"
        samtools merge $mergedName ${mergedsrrNames}
        wait $!
        samtools view -h -o $mergedSamName $mergedName
        wait $!
        samtools sort $mergedName -o $sortedMergedName
        wait $!
        samtools index $sortedMergedName
        wait $!
        ./reading $mergedSamName $resultMergedRal
        wait $!
        python search_transposeable.py ${a[0]}
        rm $mergedSamName

    else # go here for every other normal SRR entry
        ralName=${a[1]}".bam"
        ralSamName=${a[1]}".sam"
        sortedralName=${a[1]}"_sorted.bam"
        resultral=${a[0]}"-result"
        if [ ${#a[0]} -ne 1 ]
        then
            samtools view -h -o $ralSamName ${ralName}
            wait $!
            samtools sort ${ralName} -o $sortedralName
            wait $!
            ./reading $ralSamName $resultral
            wait $!
            python search_transposeable.py ${a[0]}
            rm $ralSamName
        fi
    fi
done < "$1"
