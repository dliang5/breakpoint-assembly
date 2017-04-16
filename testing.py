"""
This is a test file to see if the TE elements are being read correctly with the matches
"""
class TE(): 
    def __init__(self, chrom, pos1, pos2) : 
        self.chrom = chrom 
        self.pos1 = pos1 
        self.pos2 = pos2
    
    def getChrom(self): 
        return self.chrom
    def getPos1(self): 
        return self.pos1 
    def getPos2(self): 
        return self.pos2

class Cluster(): 
    def __init__(self, chrom, lpos1, lpos2, hpos1, hpos2): 
        self.chrom = chrom 
        self.lpos1 = lpos1 
        self.lpos2 = lpos2 
        self.hpos1 = hpos1 
        self.hpos2 = hpos2 

    def getChrom(self): 
        return self.chrom 
    def getlPos1(self): 
        return self.lpos1
    def getlPos2(self): 
        return self.lpos2
    def gethPos1(self): 
        return self.hpos1 
    def gethPos2(self): 
        return self.hpos2
TE_store = [] 
clu_Store = [] 

with open("testingTE.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split("\t")]
        decoy = TE(content[0], int(content[3]), int(content[4]))
        TE_store.append(decoy)
with open("summary-set_5_8-result", 'r') as f: 
    for i in f: 
        if i in ['\n', '\r\n']:
            continue 
        content = [x.strip() for x in i.split("\t")]
        decoy = Cluster(content[2], int(content[3]), int(content[5]), int(content[4]), int(content[6]) )
        clu_Store.append(decoy) 

# checking to see if they match or not and write into the new file.
        
