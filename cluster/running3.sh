# this is the third version of running2.sh which uses officialScript.py 
#!/bin/bash 

# this should do all of the processes from downloading the SRR file to running the assembly.
srrName=$2
srrName1=$3 # this is optional variable
setName=$1
setFile="${setName}.bam"
name1Fastq1="${srrName1}_1.fastq.gz" 
name1Fastq2="${srrName1}_2.fastq.gz" 
name1BamName="${srrName1}.bam"

fastq1="${srrName}_1.fastq.gz" 
fastq2="${srrName}_2.fastq.gz"
bamName="${srrName}.bam" 

sortedBam="${setName}_sorted.bam"
samName="${setName}.sam"


resultFile="${setName}-result"
# breakPointFile="${setName}-breakpoints"

fastq-dump --gzip --split-files $srrName
bwa mem ../../DGRP/1_SRR_file/chrom_ref $fastq1 $fastq2 | samtools view -Sb - > $bamName

# checking for a second SRRname file here to append
if [ $srrName1 -eq 0 ]
then 
    samtools view -h -o $samName $bamName
    samtools sort $bamName -o $sortedBam 
    samtools index $sortedBam 

    #done downloading and mapping everything here now 
    # idBreakPoints="id${breakPointsFile}.txt" #getting the ID names 
    # BreakPointReads="${setName}_reads.txt" #obtained all the reads from fastq 
    ./reading $samName $resultFile
    python search_transposeable.py $setName  
else
    # second SRRname is found 
    fastq-dump --gzip --split-files $srrName1 
    bwa mem ../../DGRP/1_SRR_file/chrom_ref $name1Fastq1 $name1Fastq2 | samtools view -Sb - > $name1BamName
    samtools merge $setFile $bamName $name1BamName
    samtools view -h -o $samName $setFile 
    samtools sort $setFile -o $sortedBam
    samtools index $sortedBam 
    ./reading $samName $resultFile 
    python search_transposeable.py $setName
fi 

#cleaning up here aka removing all the sam and fastq files 
rm $name1Fastq1 $name1Fastq2 $setName $fastq1 $fastq2 $samName

# if [-s $breakPointFile ]
# then 
#     echo "${srrName} has nothing in it" >> badBreakPoints
#     #badBreakPoints = nothing of potential breakpoints 
# # else 
# #     asFile="${srrName}_1.out" #file produced from fastq2fa_qual.pl in fasta format 
# #     qualAsFile="${srrName}_1.out.qual" #file produced from fastq2fa_qual.pl in qual mode 
# #     cut -f2 < $breakPointsFile > $idBreakPoints
# #     gzip -dc $fastq1 | perl parse_reads.pl $idBreakPoints > $BreakPointReads
# #     gzip -dc $fastq2 | perl parse_reads.pl $idBreakPoints >> $BreakPointReads
# #     perl fastq2fa_qual.pl < $BreakPointReads > $asFile 2> $qualAsFile 
# #     # phrap -vector_bound 0 -forcelevel 10 $BreakPointReads
# #     echo "${srrName} is done"
# fi 
    