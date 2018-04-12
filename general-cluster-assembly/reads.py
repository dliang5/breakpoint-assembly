""" this will search through the fastq files and get the selected IDs and what not """ 

import os , sys 

currentID = '' 
idList = []
with open(sys.argv[1], 'r') as f: 
    for line in f: 
        idList.append(line.strip("\n"))

setflag = 0
currentId = ""
# this checks for the id first if there's an id it'll set a flag 
# this will then automatically set the next line in right away
for index, line in enumerate(sys.stdin):
    #print(index, index % 4, line) 
    content = line.strip("\n")
    
    if '@' in content and 'SRR' in content: 
        curid = content.split(" ")[0][1:]
        if curid in idList: 
            currentId = curid
            setflag = 1 
            print(content)
    elif setflag == 1 and index % 4 == 1 and curid == currentId: 
        print(content)
    elif setflag == 1 and index % 4 == 2 and curid == currentId: 
        print(content)
    elif setflag == 1 and index % 4 == 3 and curid == currentId: 
        print(content) 
        setflag = 0
