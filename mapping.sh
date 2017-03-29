#!/bin/bash 
#$ -cwd
#$ -j y 
#$ -S /bin/bash 
# the idea is to do the entirity of the SRR mapping in one go. 
COUNTER=1
SET=$2
while IFS='' read -r line || [[ -n "$line" ]]; do 
    line1="${line}_1" 
    line2="${line}_2"   
    name1="set_${SET}_${COUNTER}.bam" 
    name2="set_${SET}_${COUNTER}_sorted.bam"
    name3="set_${SET}_${COUNTER}.sam"
    fastq1="${line1}.fastq.gz"
    fastq2="${line2}.fastq.gz"
    fastq-dump --gzip --split-files $line 
    wait $!
    bwa mem ../../DGRP/1_SRR_file/chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $name1 
    wait $!
    samtools view -h -o $name3 $name1 
    wait $!
    samtools sort $name1 -o $name2 
    wait $!
    samtools index $name2
    wait $!
    let COUNTER=COUNTER+1  
done < "$1" 
# fastq-dump --gzip --split-files $line
# bwa mem ../../DGRP/1_SRR_file/chrom_ref SRR654649_1.fastq.gz SRR654649_2.fastq.gz | samtools view -Sb - > set_13_1.bam
# samtools view -h -o set_13_1.sam set_13_1.bam
# samtools sort set_13_1.bam -o set_13_1_sorted.bam
# samtools index set_13_1_sorted.bam
