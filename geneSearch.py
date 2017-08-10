
# a simple gene search
import sys
from random import randint 
# requires two entries
class inversion:
    def __init__(self, chromRef,lpos1, lpos2, hpos1, hpos2, offset, highOffset):
        self.chromRef = chromRef
        self.lpos1 = int(lpos1)
        self.lpos2 = int(lpos2)
        self.hpos1 = int(hpos1)
        self.hpos2 = int(hpos2)
        self.offset = int(offset)
        self.highOffset = int(highOffset)

    def getlpos1(self): 
        return (self.lpos1)
    def getlpos2(self): 
        return (self.lpos2) 
    def gethpos1(self): 
        return (self.hpos1)
    def gethpos2(self): 
        return (self.hpos2) 
    def getoffset(self):
        return (self.offset)
    def gethighoffset(self): 
        return (self.highOffset) 
    def getChromRef(self): 
        return (self.chromRef) 
    def display(self): 
        print ("{} {} {} {} {} {} {}".format(self.chromRef, self.lpos1, self.lpos2, self.hpos1, self.hpos2, self.offset, self.highOffset))

# used for dmel-gene.txt 
class gene:
    def __init__(self, chromRef, pos1, pos2):
        self.chromRef = chromRef
        self.pos1 = int(pos1)
        self.pos2 = int(pos2)

    def getChromRef(self):
        return (self.chromRef) 
    def getpos1(self): 
        return (self.pos1) 
    def getpos2(self): 
        return (self.pos2)
    def display(self): 
        print("{} {} {}".format(self.chromRef, self.pos1, self.pos2)) 

# just a placeholder here  
class pos: 
    def __init__(self, pos, gene=0):
        self.pos = pos 
        self.gene = gene

    def getpos(self):
        return int(self.pos) 
    def setGene(self, genePos): # might have a list for two entries to compare against. 
        self.gene = genePos 
    def getgene(self): 
        return int(self.gene) 


# this checks the gene list and get the first position of it to compare against 
# another first position. Might compare the second position as well in case tho. 
def checkGeneList(selectedPos, gene): 
    # should i get the chrom arm based? 
    if gene.getpos1() <= selectedPos.getpos() <= gene.getpos2():
        selectedPos.setGene(gene.getpos1()) # keeping an easy number to compare again 

# starting point here for the rest of the script ---------------------
chromList = ["3L", "2R", "3R", "2L", "X"]
chromDict = {"3R": 32079331, "3L": 28110227, "2R":25286936, 
             "X": 23543271, "2L": 23513712}
inversionList = list() # storing all imports here
geneList = list()
hit = 0 
count = 0 
# parsing the inversion file out
with open("Breakpoints.csv", 'r') as inversionEntry:
    while True:
        entry1 = inversionEntry.readline()
        if(entry1.split(",")[0] == "Inversion "): continue # getting rid of the Inversion line and restarting
        entry2 = inversionEntry.readline()
        if not entry2 or "in betweens" == entry1.split(",")[0] == "in betweens": break

        content1 = entry1.split(",")
        content2 = entry2.split(",")
        if content1[0] == "In between": # not going to the common inversion, stopping here
            break 
        chrom = (content1[0].split("("))[1].split(")")[0]
        # print (content1)
        # print (content2)
        # print(chrom) 
        invClass = inversion(chrom, content1[1], content1[2], content2[1], content2[2], content1[3], content2[3] )
        inversionList.append(invClass)

# parsing the gene file out here
with open("dmel-gene.txt", 'r') as geneEntry:
    for line in geneEntry:
        content = line.split("\t")
        if(content[0] not in chromList): continue
        geneClass = gene(content[0], content[3], content[4])
        geneList.append(geneClass)
        # geneClass.display() # testing to see if everything went through as planned.   
for i in range(1000):
    hit = 0 
    count = 0 
    already_list = set()
    # this is where the permutation test begins.
    for invEntry in inversionList: 
        low1 = pos(randint(0, chromDict[invEntry.getChromRef()] - invEntry.gethpos2()))
        low2 = pos(low1.getpos() + invEntry.getoffset() )
        high1 = pos(low1.getpos() + invEntry.gethighoffset())
        high2 = pos(low1.getpos() + invEntry.gethpos2() )
        
        # match the gene and get the first position of said gene
        # gotta check the entire gene list and keep on adding it up
        for gene in geneList: 
            checkGeneList(low1, gene) 
            checkGeneList(low2, gene) 
            checkGeneList(high1, gene) 
            checkGeneList(high2, gene) 

            if low1.getgene() in already_list or low2.getgene() in already_list or high1.getgene() in already_list or high2.getgene() in already_list:
                continue
            # comparing to see if they have matching gene positions or not
            # print("{} vs {}".format(low1.getgene(), low2.getgene()))
            # print("{} vs {}".format(high1.getgene(), high2.getgene()))
            if(low1.getgene() == low2.getgene()): 
                already_list.add(low1.getgene())
                hit += 1 
            if (high1.getgene() == high2.getgene()): 
                already_list.add(high1.getgene())
                hit += 1
            count += 2
    print("hit count is {}".format(hit))