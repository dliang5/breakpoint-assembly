## Everything outside of the folders are the tools needed to cluster. 
1. cluster.py (and reading3_1.cpp once it's fixed) are the files that cluster reads
    
    1. virutually does the same thing (reading3_1.cpp and cluster.py)

2. search_trans.py - post processing for dmel files and matches TE with the breakpoints to ensure more positive inversions

3. fastq2fa_qual.pl and parse_reads.pl - converts fastq to fasta and quality and parse the ID. 

4. running3-1.sh - is how I ran the entire program save for phrap assembly to not waste space. 

## How some programs are ran for individual purposes: 

1. ```python2.7 or python3 cluster.py < *.sam > (ex. python3 cluster.py < 857.sam > )```
   output : good_<*>-result - full view of the clusters with both forward and reverse clusters 
            summary_<*>-result 

2. ```python2.7 search_trans.py * (ex. python2.7 search_trans.py 857)``` 

3. ```perl fastq2fa_qual.pl < $fastq > $1 2> $2``` 
    1.  ex. ```perl fastq2fa_qual.pl < ${fastq_file} > ${fasta_file} 2> ${quality_file}```

4. ```parse_reads.pl```  - the input are piped using gzip -dc of the fastq file 
    1.  use "cut -f2 <target-id-set>" 
    2.  ex. gzip -dc $SRRfastq | perl parse_reads.pl $idSUBMITFILE > or >> $fastq_file 

5. ```phrap - phrap -vector_bound 0 -forcelevel 10 $fasta_name (quality file has to be of the same name)``` 



## 3rd Party tools and Assumptions: 
1. fastq-dump 
2. bwa mem
3. phrap - phrap -vector_bound 0 -forcelevel 10 $fasta_name (quality file has to be of the same name) 
