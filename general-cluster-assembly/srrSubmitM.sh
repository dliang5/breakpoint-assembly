#!/bin/bash
#SBATCH --job-name -test-job
#SBATCH --nodes 1
#SBATCH --output output.out
fq1=$1"_1.fastq.gz"
fq2=$1"_2.fastq.gz"
bamFile=$1".bam"
samFile=$1".sam"
fastq-dump --gzip --split-files $1
wait $!
# chrom_ref = the index of the chromosome regions
bwa mem chrom_ref $fq1 $fq2 | samtools view -Sb - > $bamFile
wait $!

fq11=$2"_1.fastq.gz"
fq22=$2"_2.fastq.gz"
bamFile1=$2".bam"
samFile1=$2".sam"
fastq-dump --gzip --split-files $2
wait $!
# chrom_ref = the index of the chromosome regions
bwa mem chrom_ref $fq11 $fq22 | samtools view -Sb - > $bamFile1
wait $!
samtools merge $3".bam" $bamFile $bamFile1

samtools view -h -o $3".sam" $3".bam"
wait $!

### clustering here 
python2.7 cluster.py $3".sam"
wait $!
rm $3".sam" # to save disk space
