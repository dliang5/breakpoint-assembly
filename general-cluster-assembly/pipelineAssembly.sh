#!/bin/bash
#SBATCH --job-name -mosquito-job
#SBATCH --nodes 1 
#SBATCH --output mosquitotest.out
#SBATCH --ntasks=3

# this calls the other five files to read and parse the files for assembly 
# after this then we are done for now 


# what do I have to work on this

#1. call parse_id.py to obtain the id from the bam file
#2. call reads.py to gather the fastq files
#3. then call fastq2fa_qual.pl to obtain the fasta and qual files 
#4. then and finally call the assembly of both of them and you're done

bam=$1
chrom=$2
firstpos=$3 
secondpos=$4
#1 and #2 here 
samtools view $bam $chrom":"$firstpos"-"$secondpos | python2.7 parse_id.py $chrom $firstpos $secondpos > $1".id" 
#3 
perl fastq2fa_qual.pl < $1".id" > $1".fasta" 2> $1".fasta.qual"
#4
phrap -vector_bound 0 -forcelevel 10 $1".fasta"