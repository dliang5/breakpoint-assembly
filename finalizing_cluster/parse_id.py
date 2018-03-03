""" this is a file that produces the id of the region inputted from samtools -r 
    this should work with the perl file to produce a more defined breakpoint compared to 
    doing it manually

    steps to improve: let user add inputs to run it and find those particular breakpoints

    how this works currently: samtools view <BAM_file> -r REGION | python2.7 parse_id.py
    ex. DGRP-857 has a inversion breakpoint pair at 18mil and 12mil
    samtools view SRR835063.bam -r chr2R:18000000-19000000 | python2.7 parse_id.py
        <<< this does not work: new works - 
    samtools view SRR835063.bam | python2.7 parse_id.py <chrom> <first_location> <second_location>
    ex. samtools view SRR835063.bam | python2.7 parse_id.py 2R 18836300 12515000
"""


import sys, os
# taking the input of the samtool region and put them into the file
outputfile = ''
option = ''
a = sys.stdin.readline() 
checkingID = (a.split("\t")[0]).split(".")[0] 
isFileThere = checkingID+"_id.txt" 

chrLoc = sys.argv[1]
firstLoc = int(sys.argv[2])
secondLoc = int(sys.argv[3]) 

if os.path.exists(isFileThere): 
    option = 'a'
else: 
    option = 'w' 

# adding the first id line back into there
with open(isFileThere, option) as wf:
    for line in sys.stdin: 
        content = line.split("\t")
        chr = content[2] 
        first = int(content[3])
        second = int(content[7])

        if chr == chrLoc: 
            # check if they are within range
            if (firstLoc - 1000) < first and (firstLoc + 1000) > first and (firstLoc - 1000) < second and (firstLoc + 1000) > second: 
                wf.write(content[0] + "\n") 
            if (secondLoc - 1000) < second and (secondLoc + 1000) > second and (secondLoc - 1000) < first and (secondLoc + 1000) > first: 
                wf.write(content[0] + "\n") 
