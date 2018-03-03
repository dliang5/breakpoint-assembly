""" this is an easier way to search for the corresponding genome to its SRR counterpart 
    INPUT: python2.7 name
    OUTPUT: <name's SRRname> 
"""

import sys, subprocess as sp


entry = sys.argv[1]
# note: Ral-857's SRX021492 works here now check for ral-857 as well.
entry = entry.upper()
newEntry = ""
found = False
# checks if it's actually RAL-857
# if is, go bring it from "TableS1_individuals.csv"
if "SRX" not in entry or "ERX" not in entry: 
    with open("TableS1_individuals.csv", 'r') as tableFile:
        for tableEntry in tableFile:
            content = tableEntry.split(",")
            if entry in content[0]:
                # print(content)
                newEntry = content[3]
                if ("\"" in newEntry):
                    newEntry = content[4].strip("\"")
                found = True
                break 

if found == True:
    entry = newEntry
    if "SRR" in entry or "ERR" in entry:
        print(entry)
    else:
        cmd = "grep -r " + entry.rstrip("\n") + " SRA_Accessions.tab"
        content = sp.Popen(
            cmd, stdout=sp.PIPE, shell=True
        )
        out = str(content.communicate()[0])
        check1 = out.split("\t")[0]
        if "SRR" in check1 or "ERR" in check1: 
            # return check1 back, otherwise return 0 or something
            print(check1)
        else: 
            print("NONE") # this will make align.sh skip/continue to next iteration
else:
    print("NONE")
