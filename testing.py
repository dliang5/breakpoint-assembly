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
    def __init__(self, index, chrom, lpos1, lpos2, hpos1, hpos2, size): 
        self.index = index 
        self.chrom = chrom 
        self.lpos1 = lpos1 
        self.lpos2 = lpos2 
        self.hpos1 = hpos1 
        self.hpos2 = hpos2 
        self.size = size 
    def getIndex(self): 
        return self.index
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
    def getSize(self): 
        return self.size
TE_store = [] 
clu_Store = [] 

with open("TE_sites.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split("\t")]
        decoy = TE(content[0], int(content[3]), int(content[4]))
        TE_store.append(decoy)
        
with open("summary-set_5_8-result", 'r') as f: 
    for i in f: 
        if i in ['\n', '\r\n']:
            continue 
        content = [x.strip() for x in i.split("\t")]
        decoy = Cluster(content[0], content[2], int(content[3]), int(content[5]), int(content[4]), int(content[6]), int(content[7]) )
        clu_Store.append(decoy) 

print (len(clu_Store)) 
print (len(TE_store ))
tranClus = []
# checking to see if they match or not and write into the new file.
for i in clu_Store: 
    counter = 0
    for j in TE_store: 
        total = (j.getPos1() + j.getPos2()) / 2 
        if i.getChrom() == j.getChrom(): 
            if (j.getPos1() <= i.getlPos1() <= j.getPos2()): 
                tranClus.append(i) 
                print("{0} <= {1} <= {2} hi".format(j.getPos1(), i.getlPos1(), j.getPos2())) 
                break 
            elif (j.getPos1() <= i.getlPos2() <= j.getPos2()): 
                tranClus.append(i)  
                print("{0} <= {1} <= {2} hi1".format(j.getPos1(), i.getlPos2(), j.getPos2())) 
                break 
            elif (j.getPos1() <= i.gethPos1() <= j.getPos2()): 
                tranClus.append(i)   
                print("{0} <= {1} <= {2} hi2".format(j.getPos1(), i.gethPos1(), j.getPos2()))               
                break 
            elif (j.getPos1() <= i.gethPos2() <= j.getPos2()): 
                tranClus.append(i)      
                print("{0} <= {1} <= {2} hi3".format(j.getPos1(), i.gethPos2(), j.getPos2()))            
                break
            # trying to find a cluster within the range
            else:
                pass 
                # print("{0} vs {1}".format(counter, i.getSize()))

for i in tranClus: 
    Status = False 
    for j in clu_Store: 
        if i.getIndex() == j.getIndex(): 
            Status = True # there exist an element of that index 
    clu_Store.remove(i)

# write_file = "set_5_8-test-file.txt" 
# writing = open(write_file, 'w') 
# for i in clu_Store: 
#     with open("good_ZI213-result", 'r') as f: 
#         for line in f: 
#             content = [x.strip() for x in line.split("\t")]
#             if content[0] == i.getIndex(): 
#                 # then start writing everything here
#                 writing.write(line) 


