import os, sys, sp as subprocess
inputfile = sys.argv[1] # the ERS and com file 
comp = sys.argv[2] # SRA_accession.tab
with open(inputfile, 'r') as f: 
    for entry in f: 
        content = (entry.strip("\n")).split(",")
        ers = content[1] 
        # right now, use sp for grep search and pick the last ERR file of that search 
        # then use that to collect the others as well
