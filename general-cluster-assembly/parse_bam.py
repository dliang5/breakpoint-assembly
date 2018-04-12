""" 
    This program will take in the lines from a BAM file from the command, and it 
    is supposed to function similar to samtools view region -r if it ever comes 
    down to the case it can't be used for whatever reasons.

    input: output of bam to sam file
        to run: samtools view -h <SRR.bam> | python2.7 parse_bam.py <chr> <breakpoint1> <breakpoint2>  
    output: SRRid of a region about
"""

import sys, os
chrIn = sys.argv[1] 
breakIn1 = int(sys.argv[2])
breakIn2 = int(sys.argv[3]) 

for line in sys.stdin: 
    content = (line.strip("\n")).split("\t") 
    if len(content) < 7:
        continue
    chr = content[2]
    break1 = int(content[3]) 
    break2 = int(content[7]) 
    if chr != chrIn: 
        continue 
    # gonna use a 1kb window size to make buffer and go around and see if there are other breakpoints.
    # might try out 100kb window size later on 
    elif (breakIn1 - 100000) < break1 and (breakIn1 + 100000) > break1 and \
        (breakIn2 - 100000) < break2 and (breakIn2 + 100000) > break2: 
            print(content[0])
