#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
# qsub align.sh overall-1-file
# this is just to download and align files that's all up to bam file.
while IFS='' read -r line|| [[ -n "$line" ]]; do
    # echo "text from file " $line
    a=( $line )
    # {a[0]} = name of file and {a[1...n]} = SRRnames and erx ones
    srrNames1="${a[0]}"
    srrNames=""
    mergedsrrNames=""
    # for word in $line; do echo $word; done
    # echo "here are the individual stuff in it" "${a[1]}"

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
            #fastq-dump --gzip --split-files $word
            wait $!
            #bwa mem ../../../../reference/chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $bamName
            wait $!
        fi
    done
    b=( $srrNames )
    # if [ "${#b[@]}" -gt 1  ]
    # then
    #     fastq1="${myArgs[i]}_1.fastq.gz"
    #     fastq2="${myArgs[i]}_2.fastq.gz"
    #     bamName="${myArgs[i]}.bam"
    #     sortedBam="${myArgs[i]}_sorted.bam"
    #     fastq-dump --gzip --split-files "${myArgs[i]}"
    # else
    #     echo "NOPE"
    # fi

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
    else
        #for i in "${!myArgs[@]}"; do
        # variables to be used or something
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
        #done
    fi

done < "$1"
