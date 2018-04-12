""" this will output a file that contaians all of the lines that does match the ones in the csv files """ 
import os, sys
comp1 = sys.argv[1] # this is the file im given to compare 
comp2 = sys.argv[2] # this is the entire population list 

with open(comp1, 'r') as f: 
    for entry in f: 
        content = (entry.strip("\n")).split(",") 
        name = content[5] 
        # then comparing to the other file here
        with open(comp2, 'r') as cpf: 
            for compLine in cpf: 
                litLine = (compLine.strip("\n")).split()
                litname = content[1]
                if name == litname or name in litLine: 
                    print("{},{}".format(litLine[1], litLine[2]))
                    continue 
