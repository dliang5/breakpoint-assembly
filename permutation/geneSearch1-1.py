""" Permutation test for rare inversion breakpoints with common inversion breakpoints as the null hypothesis
The main idea is that rare inversion should break or "appear" more often than common inversion breakpoints as
they are deletroious based on observation.

This program reads both rare and common inversion breakpoints data and the full gene expression list.
For each breakpoint, it will randomly select a position within its respective chromosome arm
then it will compare if there was an actual match at that location on the gene list up to 1000+ iterations.
If both forward and reverse match (both low and high positions) then it will count.

Expectation: the percentage is lower  on rare than common
             or - calculating the means for all of the iterations and getting the differences
             if it is less than 5% p-value than we are good to go.
"""
# a simple gene search
import sys
from random import randint
# requires two entries
class inversion:
    def __init__(self, chromRef,lpos1, lpos2, hpos1, hpos2, offset):
        self.chromRef = chromRef
        self.lpos1 = int(lpos1)
        self.lpos2 = int(lpos2)
        self.hpos1 = int(hpos1)
        self.hpos2 = int(hpos2)
        self.offset = int(offset)
        # self.highOffset = int(highOffset)

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
    # def gethighoffset(self):
    #     return (self.highOffset)
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
def checkGeneList(selectedPos, gene, chrom):
    if gene.getChromRef() == chrom and gene.getpos1() < selectedPos.getpos() < gene.getpos2():
            #print("-----------selecting gene {}------------".format(gene.getpos1()))
            selectedPos.setGene(gene.getpos1()) # keeping an easy number to compare again
    else: 
        selectedPos.setGene(0)

# starting point here for the rest of the script ---------------------
chromList = ["3L", "2R", "3R", "2L", "X"]
chromDict = {"3R": 32079331, "3L": 28110227, "2R":25286936,
             "X": 23543271, "2L": 23513712}
inversionList = list() # storing all imports here
geneList = list()
hit = 0
count = 0

# parsing the inversion file out for both common and rare inversions
with open("Breakpoints3.csv", 'r') as inversionEntry:
    # getting rare inversion breakpoints here
    for times in range(2):
        currentList = list() # getting the first list for rare inversion
        while True:
            entry1 = inversionEntry.readline()
            if entry1.split(",")[0] == "Inversion ": continue # getting rid of the Inversion line and restarting
            if "In between" in entry1: break # moving on to the common inversions
            entry2 = inversionEntry.readline()
            content1 = entry1.split(",")
            content2 = entry2.split(",")
            # this if stmt is prob bad since the one above already checks for it
            if content1[0] == "In between" or content1[0] == "Unique set": break
            chrom = (content1[0].split("("))[1].split(")")[0]
            #print(content2[3])
            invClass = inversion(chrom, content1[1], content1[2], content2[1], content2[2], content1[3])
            # print("{} {} {} {} {} {}".format(chrom, content1[1], content1[2], content2[1], content2[2], content1[3]))
            currentList.append(invClass)
        inversionList.append(currentList)

# parsing the gene file out here
with open("dmel-gene.txt", 'r') as geneEntry:
    for line in geneEntry:
        content = line.split()
        if(content[0] not in chromList): continue
        name = content[8].split(";")
        if "CR" in name[1]: continue
        geneClass = gene(content[0], content[3], content[4])
        geneList.append(geneClass)
        # geneClass.display() # testing to see if everything went through as planned.

# this is to store each iteration in their respective spot for mean computation purposes
meanList = dict()
graphList = dict()

for i in range(2):
    meanList[i] = 0
    graphList[i] = list()
    
# permutation test here. Iterating at the min. 1000 times.
for j in range(2):
    targetCount = 0
    if j == 0:
        targetCount = 20 # rare inversion count score
    elif j == 1:
        targetCount = 4 # common inversion count score
    # remove everything up from this area to keep it back the original code style.
    # starting permutation iteration here

    """ 
    try modifying the script to count the breakpoints for each
    """
    for i in range(1000):
        hit = 0
        # this is where the permutation test begins.
        for index, invEntry in enumerate(inversionList[j]):
            low1 = pos(randint(0, chromDict[invEntry.getChromRef()] - abs(invEntry.gethpos2() - invEntry.getlpos1() )))
            low2 = pos(low1.getpos() + invEntry.getoffset() )
            high1 = pos(low1.getpos() + abs(invEntry.gethpos1() - invEntry.getlpos1() ))
            high2 = pos(low1.getpos() + abs(invEntry.gethpos2() - invEntry.getlpos1() ))

            hitlow = 0
            hithigh = 0
            # match the gene and get the first position of said gene
            # gotta check the entire gene list and keep on adding it up
            for gene in geneList:
                checkGeneList(low1, gene, invEntry.getChromRef())
                checkGeneList(low2, gene, invEntry.getChromRef())
                checkGeneList(high1, gene, invEntry.getChromRef())
                checkGeneList(high2, gene, invEntry.getChromRef())
                # comparing to see if they have matching gene positions or not
                if(low1.getgene() == low2.getgene() and low1.getgene() != 0 and low2.getgene() != 0):
                        hitlow += 1

                if (high1.getgene() == high2.getgene() and high1.getgene() !=0 and high2.getgene() != 0):
                        hithigh += 1
            if ( hitlow > 0 or hithigh > 0 ):
                hit+= 1
            # if ( hithigh > 0):
            #     hit += 1
        # graphing purpose
        graphList[j].append(hit)
        # for the p-value
        if hit > targetCount:
            meanList[j] += 1
# getting the results here
for key, value in meanList.items():
    print("{} -> {} {}".format(key, value, ((value * 100)/1000)*100))
    print("{} -> {} {}".format(key, value, value/float(1000)))

for i in range(2):
    for j in graphList[i]:
        print("{} {}".format(i, j))
