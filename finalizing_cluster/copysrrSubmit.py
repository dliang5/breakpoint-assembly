""" this is an automation for srrSubmit.sh
    input: python2.7 srrSubmit.py <srrText>
    output: good_<srr>-breakpoints summary_<srr>-breakpoints """


""" description: 
    This runs on known files and pipeline the shit in 
    
    In progress: 
    Create a portion that checks if args is 2 or 1
    if 1, then check the list of files and keep on running the job nonstop 
    else 2, then compare and find the look. 

    """ 
import os, sys, subprocess as sp

# this check if there are SRX or multiple SRX.. vice versa with SRR
def check(content): 
    tempList = list()
    flag = 0
    for i in content: 
        if "ERX" in i or "ERR" in i:
            flag = 1
    for ele in content: 
        if "SRX" in ele or "ERX" in ele: 
            # run the sb process here and wait for it to finish
            ele1 = ele 
            if "\"" in ele: 
                ele1 = ele.replace("\"", "") 
            grep = "grep -r {} SRA_Accessions.tab".format(ele1)
            result = sp.check_output(grep, shell=True)
            resultCheck = result.split("\t")
            if flag == 0:
                if "SRR" not in resultCheck[0]: 
                    continue 
            elif flag == 1: 
                if "ERR" not in resultCheck[0]: 
                    continue
            tempList.append(resultCheck[0])
        elif "SRR" in ele or "ERR" in ele:
            ele1 = ele
            if "\"" in ele: 
                ele1 = ele.replace("\"", "")
            
            tempList.append(ele1) 
    return tempList 

nameSRR = dict()
checkSRR = list()
# opening the file of known inversions or checking known inversions already found  

with open(sys.argv[2], 'r') as f: 
    for line in f: 
        if line == "": 
            continue
        content = (line.strip("\n")).split() 
        checkSRR.append(content[0])

# opening the file of john pool
with open(sys.argv[1], 'r') as f: 
    for line in f: 
        content = (line.strip("\n")).split(",")
        if content[0] not in checkSRR: 
            continue 

        tempList = check(content)
        nameSRR[content[0]] = tempList
print(nameSRR) 

for k,v in nameSRR.items(): 
    if len(v) == 1: 
        run = "sbatch srrSubmit.sh {} {}".format(v[0], k) 
        # run the script here 
        sp.Popen(run, shell=True, exectuable='/bin/bash')
    else: 
        run = "sbatch srrSubmitM.sh"
        for i in v: 
            run += " {}".format(i)
        run += " {}".format(k)
        # run the merged script here
        sp.Popen(run, shell=True, exectuable='/bin/bash')
