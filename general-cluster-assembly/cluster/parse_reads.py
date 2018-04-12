""" 
    input gzip -dc <fastq_file> | python2.7 parse_reads.py <id_file> 
    This file will parse the corresponding fastq lines and compare against <id_file>
    to see if any of them matches or not.  
"""
import sys
ID_list = list() 
with open(sys.argv[1], 'r') as id_file:
    for ID in id_file: 
        ID_list.append(ID.split(".")[1]) 
        
        
counter = 0
for line in sys.stdin:
    entry = line.rstrip("\n")
    
    if counter == 1: 
        print line
        counter = 0  
        
    if "@" in line:
        for ID in ID_list: 
            if ID in line:
                print line
                counter = counter + 1
                break
