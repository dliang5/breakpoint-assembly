#!/bin/bash
#SBATCH --job-name -test-job
#SBATCH --nodes 1
#SBATCH --output output1.out
fq1=$1"_1.fastq.gz"
fq2=$1"_2.fastq.gz"
bamFile=$2".bam"
samFile=$2".sam"
fastq-dump --gzip --split-files $1
wait $!
# chrom_ref = the index of the chromosome regions
bwa mem chrom_ref $fq1 $fq2 | samtools view -Sb - > $bamFile
wait $!
samtools view -h -o $samFile $bamFile
wait $!

### clustering here 
python2.7 cluster.py $samFile
wait $!
