## ASSUMPTION AND KNOWS


## I only know: 

1. the filename the user will input and if it is not ERR or SRR then I have to correspond that to it then
   SRA_accession.tab then

2. srrsubmit.py cannot work for this case and would need a redesign or another creation of it
    1. srrsubmit.py reads the names and searches the corresponding
        as such. So we would have to modify it to include poplulation.csv. 
    2. once that is done, create a pipeline program (a general one)        that could handle all types please this time.
    3. I think maybe you should consider that before working on            mosquito

3. The user will run 
    ```sbatch testing.sh {files to look at and cluster}```

4. The assumption is that they have a SRA_accession.tab nearby to finalize and double check the results 
    1. TODO: (backlog) - retrieve the data from SRA online instead to save space -> that would be for the other cases

5. So when converting between ERS and SRS convert to the last one as        that is generally the bigger file so perhaps better accuracy? no       lie
