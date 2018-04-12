""" this is a python file that only parse the id of the selection breakpoints and their first and second
locations just a much more smaller program than parse_id.py"""


import sys, os 

userChr = sys.argv[1] 
firstLoca = int(sys.argv[2]) 
secondLoca = int(sys.argv[3]) 
idList = []


for line in sys.stdin:
    line1 = line.strip("\n")
    if line1 == "": 
        continue 
    content = line1.split("\t") 
    first = int(content[5])
    second = int(content[6]) 
    chr = content[3] 
    if chr == userChr: 
        # checking to see if the first position is within the window size 
        if (firstLoca - 1000) < first and (firstLoca + 1000) > first: 
            if (secondLoca - 1000) < second and (secondLoca + 1000) > second: 
                idList.append(content[1])

for i in idList: 
    print(i) 
