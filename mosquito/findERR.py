# should probably combine this with the cluster fly

import os, sp as subprocess, sys
mosqFile = sys.argv[1] 
with open(mosqFile, 'r') as f: 
    for entry in f: 
        # start searching for the corresponding ERR and ERS from the COM files
        if entry == "" or entry == None: 
            break # break out in case there is nothing at all
        content = entry.strip("\n") 