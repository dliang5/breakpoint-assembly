#!/bin/bash
#SBATCH --job-name -mosquito-job
#SBATCH --nodes 1 
#SBATCH --output mosquitotest.out
#SBATCH --ntasks=3

module load samtools/

srr=$1 #srr name of the file so i can call fastq-dump ERR or SRR not ERS 
expName=$2 #the name of the experiment src_code here 
fq1=$2"_1.fastq.gz"
fq2=$2"_2.fastq.gz" 
bamFile=$1"-"$2".bam"
samfile=$2"-"$2".sam"


#### make sure you got the mosquito chrom_ref down here before going any further 
fastq-dump --gzip --split-files $2
bwa mem chrom_ref $fq1 $fq2 | samtools view -Sb - > $bamFile
samtools view -h -o $samFile $bamFile 
python2.7 cluster.py $samFile
python2.7 search_trans.py $1"-"$2

#### done with the clustering here