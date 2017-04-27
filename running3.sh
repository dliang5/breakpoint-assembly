# this is the third version of running2.sh which uses officialScript.py 
#!/bin/bash 

# this should do all of the processes from downloading the SRR file to running the assembly.
srrName=$1
fastq1="${srrName}_1.fastq.gz" 
fastq2="${srrName}_2.fastq.gz"
bamName="${srrName}.bam" 
sortedBam="${srrName}_sorted.bam"
samName="${srrName}.sam"
fastq-dump --gzip --split-files $srrName
bwa mem ../../DGRP/1_SRR_file/chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $bamName
samtools view -h -o $samName $bamName
samtools sort $bamName -o $sortedBam 
samtools index $sortedBam 

#done downloading and mapping everything here now 
resultFile="${srrName}-result"
breakPointFile="${srrName}-breakpoints"
idBreakPoints="id${breakPointsFile}.txt" #getting the ID names 
BreakPointReads="${srrName}_reads.txt" #obtained all the reads from fastq 
./reading $samName $resultFile
python search_transposeable.py $srrName 
if [-s $breakPointFile ]
then 
    echo "${srrName} has nothing in it" >> badBreakPoints
    #badBreakPoints = nothing of potential breakpoints 
# else 
#     asFile="${srrName}_1.out" #file produced from fastq2fa_qual.pl in fasta format 
#     qualAsFile="${srrName}_1.out.qual" #file produced from fastq2fa_qual.pl in qual mode 
#     cut -f2 < $breakPointsFile > $idBreakPoints
#     gzip -dc $fastq1 | perl parse_reads.pl $idBreakPoints > $BreakPointReads
#     gzip -dc $fastq2 | perl parse_reads.pl $idBreakPoints >> $BreakPointReads
#     perl fastq2fa_qual.pl < $BreakPointReads > $asFile 2> $qualAsFile 
#     # phrap -vector_bound 0 -forcelevel 10 $BreakPointReads
#     echo "${srrName} is done"
fi 
    