""" Permutation test for rare inversion breakpoints with common inversion breakpoints as the null hypothesis
The main idea is that rare inversion should break or "appear" more often than common inversion breakpoints as
they are deletroious based on observation.

This program reads both rare and common inversion breakpoints data and the full gene expression list.
For each breakpoint, it will randomly select a position within its respective chromosome arm
then it will compare if there was an actual match at that location on the gene list up to 1000+ iterations.
If both forward and reverse match (both low and high positions) then it will count.

Expectation: the percentage is lower on rare than common
             or - calculating the means for all of the iterations and getting the differences
             if it is less than 5% p-value than we are good to go.

how to run: python2.7 geneSearch1.py 
Can't say python since python2.6 does not support {} in print("{}".format(x))
"""
# a simple gene search
import sys
from random import randint
""" This holds each inversion from Breakpoints.csv
    lpos1 = first low position, lpos2 = second low position, etc
    offset = distance between lpos1 and lpos2
"""
class inversion:
    def __init__(self, chromRef,lpos1, lpos2, hpos1, hpos2, offset):
        self.chromRef = chromRef
        self.lpos1 = int(lpos1)
        self.lpos2 = int(lpos2)
        self.hpos1 = int(hpos1)
        self.hpos2 = int(hpos2)
        self.offset = int(offset)

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
    def getChromRef(self):
        return (self.chromRef)
    def display(self):
        print ("{} {} {} {} {} {} {}".format(self.chromRef, self.lpos1, self.lpos2, self.hpos1, self.hpos2, self.offset, self.highOffset))

""" a bit reduandant here but just for readability sake """
# used for dmel-gene.txt or the gene list each permutation will go through 
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
# for low1, low2, high1, and high2. It's used to compared each breakpoint to each other.
class pos:
    def __init__(self, pos, gene=0):
        self.pos = pos
        self.gene = gene

    def getpos(self):
        return int(self.pos)
    def setGene(self, genePos): # might have a list for two entries to compare against.
        self.gene = genePos
    def getgene(self):
        return self.gene


# this checks the gene list and gets the first and second positions of each gene to compare against
# selectedPos position. If selectedPos's position overlaps between first and second of each gene 
# then it will set "setGene" to the first gene position, if that does not hit hte gene, gene is set to 0
def checkGeneList(selectedPos, gene, chrom):
    if ( gene.getChromRef() == chrom and gene.getpos1() < selectedPos.getpos() < gene.getpos2() ) :
        selectedPos.setGene(gene.getpos1()) # keeping an easy number to compare again
    else :
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
    for times in range(2):
        currentList = list() # getting the first list for the inversions
        while True:
            entry1 = inversionEntry.readline() # reading first breakpoint or the "lower" positions
            if entry1.split(",")[0] == "Inversion ": continue # getting rid of the "Inversion" line or header and restarting
            if "In between" in entry1: break # moving on to the common inversions
            entry2 = inversionEntry.readline() # reading second breakpoint or the "higher" positions
            if not entry2 or "in betweens" == entry1.split(",")[0] == "in betweens": break
            content1 = entry1.split(",")
            content2 = entry2.split(",")
            # just skips "in between" and the second 2R(mal) positoins
            if content1[0] == "In between" or content1[0] == "Unique set": break
            chrom = (content1[0].split("("))[1].split(")")[0]
            #print(content2[3])
            invClass = inversion(chrom, content1[1], content1[2], content2[1], content2[2], content1[3])
            #print("{} {} {} {} {} {}".format(chrom, content1[1], content1[2], content2[1], content2[2], content1[3]))
            currentList.append(invClass)
        inversionList.append(currentList)

# parsing the gene file out here
with open("dmel-gene.txt", 'r') as geneEntry:
    for line in geneEntry:
        content = line.split("\t")
        if(content[0] not in chromList): continue
        geneClass = gene(content[0], content[3], content[4])
        geneList.append(geneClass)

# this is to store each iteration in their respective spot for mean computation purposes
meanList = dict()
graphList = dict()

for i in range(2):
    meanList[i] = 0
    graphList[i] = list()

# permutation test here. Iterating at the min. 1000 times.
for j in range(2):
    targetCount = 0
    # [0] is for rare inversions , [1] is for common inversions
    if j == 0:
        targetCount = 23 # rare inversion count score
    elif j == 1:
        targetCount = 4 # common inversion count score

    for i in range(1):
        hit = 0
        count = 0
        # this is where the permutation test begins.
        for index, invEntry in enumerate(inversionList[j]):
        
            # has abs since the breakpoints aren't always in line from lowest to highest.
           # low1 = pos(randint(0, chromDict[invEntry.getChromRef()] - abs(invEntry.gethpos2() - invEntry.getlpos1() )))
            low1 = pos(invEntry.getlpos1())
            low2 = pos(low1.getpos() + invEntry.getoffset() )
            high1 = pos(low1.getpos() + abs(invEntry.gethpos1() - invEntry.getlpos1() ))
            high2 = pos(low1.getpos() + abs(invEntry.gethpos2() - invEntry.getlpos1() ))

            ### record whether the high adn low breakpoints were hit
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
                
                if (high1.getgene() == high2.getgene() and high1.getgene() != 0 and high2.getgene() != 0):
                        hithigh += 1
        
            #print invEntry.getChromRef(), low1.getpos(), low2.getpos(), high1.getpos(), high2.getpos(), hitlow, hithigh, hit

            if ( hitlow > 0 ) :
                hit += 1
            if ( hithigh > 0 ) :
                hit += 1
    
        # for graphing purposes
        graphList[j].append(hit)

        if hit > targetCount:
            meanList[j] += 1

# getting the results for graphing purpose 
for i in range(2):
    for j in graphList[i]:
        print("{} {}".format(i, j))
