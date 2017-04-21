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
        self.matches = False
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
    # only be used if the value of dmel applies to both areas
    def convertStatus(self): 
        self.matches = True 
    def getStatus(self): 
        return self.matches
TE_store = [] 
clu_Store = [] 
dmel_holder = [] 
# opening transposeable reference points 
with open("TE_sites.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split("\t")]
        decoy = TE(content[0], int(content[3]), int(content[4]))
        TE_store.append(decoy)
# opening the summary clusters points 
with open("summary-ZI213-result", 'r') as f: 
    for i in f: 
        if i in ['\n', '\r\n']:
            continue 
        content = [x.strip() for x in i.split("\t")]
        decoy = Cluster(content[0], content[2], int(content[3]), int(content[5]), int(content[4]), int(content[6]), int(content[7]) )
        clu_Store.append(decoy) 
# opening the cytological breakpoints table
with open("DmelMapTable.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split(",")]
        chrom_info = [x.strip() for x in content[2].split(":")]
        if (len(chrom_info) == 2): 
            pos_info = [x for x in chrom_info[1].split("..")]
            decoy = TE(chrom_info[0], int(pos_info[0] ), int (pos_info[1]) )
            dmel_holder.append(decoy)
        else: 
            continue 

print (len(clu_Store)) 
print (len(TE_store ))
print (len(dmel_holder))
tranClus = []
# checking to see if they match or not and write into the new file.
for i in clu_Store: 
    counter = 0
    for j in TE_store: 
        if i.getChrom() == j.getChrom(): 
        # How to include the range and make it correct as well? 
            if ( j.getPos1() <= i.getlPos1() <= j.getPos2() ): 
                tranClus.append(i) 
                print("{0} <= {1} <= {2} hi".format(j.getPos1(), i.getlPos1(), j.getPos2())) 
                break 
            elif ( j.getPos1() <= i.getlPos2() <= j.getPos2() ): 
                tranClus.append(i)  
                print("{0} <= {1} <= {2} hi1".format(j.getPos1(), i.getlPos2(), j.getPos2())) 
                break 
            elif ( j.getPos1() <= i.gethPos1() <= j.getPos2() ): 
                tranClus.append(i)   
                print("{0} <= {1} <= {2} hi2".format(j.getPos1(), i.gethPos1(), j.getPos2()))               
                break 
            elif ( j.getPos1() <= i.gethPos2() <= j.getPos2() ): 
                tranClus.append(i)      
                print("{0} <= {1} <= {2} hi3".format(j.getPos1(), i.gethPos2(), j.getPos2()))            
                break
            else:
                pass 
                # print("{0} vs {1}".format(counter, i.getSize()))
# removing the TE from the cluster stores
for i in tranClus: 
    Status = False 
    for j in clu_Store: 
        if i.getIndex() == j.getIndex(): 
            Status = True # there exist an element of that index 
    if Status == True: 
        clu_Store.remove(i)

final_holder = [] 
temp_deml_holder = dmel_holder  
# read in clu_Store and compare to dmel first to see if the first position match. 
for i in clu_Store: 
    for j in dmel_holder: 
        status = False
        counter = 0 
        if i.getChrom() == j.getChrom(): 
            if( j.getPos1() <= i.getlPos1() <= j.getPos2() or j.getPos1() <= i.getlPos2() <= j.getPos2()):
                # first check done, now check again on the second position to see if they actually match.
                for k in temp_deml_holder: 
                    if k.getChrom() == i.getChrom(): 
                        if(k.getPos1() <= i.gethPos1() <= k.getPos2() or k.getPos1() <= i.gethPos2() <= k.getPos2()):
                            i.convertStatus() # setting matches to true 
                            final_holder.append(i)
                            break
                if i.getStatus() == False and counter == 0: 
                     final_holder.append(i)
                     counter +=1
                     break 
# let's not focus on this and do it later on, but any single line is considered unuseable btw.
# scanning the list of final_holder
# writing everything to the summary file
with open("summary-ZI213-breakpoints.txt", 'w') as f: 
    counter = 0 
    for count, i in enumerate(final_holder): 
        if count % 2 == 0: 
            if final_holder[count+1].getChrom() == i.getChrom(): 
                if ( (final_holder[count+1].getlPos1 / i.getlPos2() > .90):  
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".\
                    format(i.getIndex(), i.getChrom(), i.getlPos1(), i.gethPos1(), i.getlPos2(), i.gethPos2(), \
                    i.getSize(), i.getStatus()) )

                    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n\n".\
                    format(final_holder[count+1].getIndex(), final_holder[count+1].getChrom(), \
                    final_holder[count+1].getlPos1(), final_holder[count+1].gethPos1(), final_holder[count+1].getlPos2(), \
                    final_holder[count+1].gethPos2(), final_holder[count+1].getSize(), final_holder[count+1].getStatus()) )
            else: 


write_file = "ZI213-test-file.txt" 
writing = open(write_file, 'w') 
for i in clu_Store: 
    with open("good_ZI213-result", 'r') as f: 
        for line in f: 
            content = [x.strip() for x in line.split("\t")]
            if content[0] == i.getIndex(): 
                # then start writing everything here
                writing.write(line) 


