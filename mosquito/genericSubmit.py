""" This is a generic submit file that is ran by testing.sh in the cluster system
    What it does: 
        1. it will read a file from testing.sh and check for either ERR and SRR
        or any other corresponding files to them in the ER or SR range
        2. Then it calls SRA_accession.tab to get the coresponding ERR and SRR 
        3. then it will run the clustering algorithm from here up to 5 times then 
        pause for 3-5 hrs (since my school is limiting me at the momemt)
""" 
import os, sys, subprocess as sp, time

identifierFile = open("identityFile.txt", 'w')
readFile = sys.argv[1] # contains all the files to be read into
with open(readFile, 'r') as f: 
    nameSRR = dict() 
    checkSRR = dict() 

    for line in f: 
        # the assumption is that the first name will be the identifier
        # find
        skipGrep = 0 #assuming that there is no ERR or SRR in the line im looking at
        content = (line.strip("\n")).split("\t")
        for i in content: 
                if "ERR" in i or "SRR" in i: 
                    skipGrep = 1
        if skipGrep != 1:
            for ele in content: 
                if "SR" in ele or "ER" in ele: 
                    if "\"" in ele: 
                        ele = ele.replace("\"", "") #removing quotes
                    grep = "grep -r {} SRA_Accessions.tab".format(ele)
                    result = sp.check_output(grep, shell=True)
                    resultCheck = result.split("\t") 
                    if "SRR" in resultCheck[0] or "ERR" in resultCheck[0]: 
                        nameSRR[content[0]] =  ele # saving the name here
        else: 
            for ele in content: 
                if "SRR" in ele or "ERR" in ele:
                    nameSRR[content[0]] = ele

    # I have all of the names and their respective SRR here as well 
    for index, k, v in enumerate(nameSRR.items()): 
        run = "sbatch" 
        writingTo = "{}".format(k) 
        if len(v) == 1: 
            run += " {} {} {}".format(k, v[0]) 
            writingTo += "\t{}".format(v)
            sp.Popen(run, shell=True) 
            identifierFile.write(writingTo)
        else: 
            run += " srrSubmitM.sh"
            for i in v: 
                run += " {}".format(i)
                writingTo += "\t{}".format(i)
            run += " {}".format(k) 
            sp.Popen(run. shell=True) 
            identifierFile.write(writingTo) 
        if (index + 1) % 5 == 0: 
            time.sleep(3600 * 3) # sleeping the program for 4 hrs to avoid usingn too much computations
writingTo.close()
