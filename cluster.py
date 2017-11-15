# !/usr/bin/python3 
# david liang (dliang5) 
import sys, copy 
""" to use it by itself - python3 Reading3.py <samfile> 
    input : samfile (ex. 857.sam)
    output : good_<sam>-result summary_<sam>-result (ex. good_857-result summary_857-result)

    this file will scan through a SAM file and create a 2D list for each chromosome up to the nth time.
    It then merges any similar clusters to avoid overfilling the output files with too much 
    results. There are 8 samflags hardcoded: 

    for1 = [65, 129, 2161, 2225] 
    rev1 = [113, 177, 2113, 2177] 

"""
""" background information and assumption:
    Can run at python 2.7 and python3.5
"""
# this is a class creation for the purpose of using the clustering algorithm 
class location(): 
    def __init__(self, id, samFlag, chromRef, pos1, mapQual, equiva, pos2): 
        self.id = id 
        self.samFlag = samFlag 
        self.chromRef = chromRef 
        self.pos1 = pos1 
        self.mapQual = mapQual 
        self.equiva = equiva 
        self.pos2 = pos2 
        self.reverse = False 
        self.forward = False 
    def getID(self): 
        return str(self.id)
    def getSamFlag(self): 
        return int(self.samFlag) 
    def getChromRef(self): 
        return (self.chromRef) 
    def getMapQual(self): 
        return (self.mapQual)
    def getPos1(self): 
        return int(self.pos1)
    def getPos2(self): 
        return int(self.pos2)
    def getEquiva(self): 
        return (self.equiva)
    def getFor(self): 
        return self.forward 
    def getRev(self): 
        return self.reverse
    def convertRev(self): 
        self.reverse = True 
    def convertFor(self): 
        self.forward = True

# this checks all of the clusters and submit the entry into that cluster within a certain range
def search(tmpList, entry, minCluster): 
    inRange = False 
    for i in tmpList: 
        if i[0].getChromRef() != entry.getChromRef(): 
            continue 
        elif entry.getPos1() < i[0].getPos1() + 3000: 
            if abs(entry.getPos2() - i[0].getPos2()) < 3000 and abs(entry.getPos1() - i[0].getPos1()) < 3000:
                i.append(entry) 
                inCluster = True
    return inRange, tmpList

# prints out the breakpoints 
def printCluster(forward, reverse, writeFile):
    for cForward in forward:
        for f in range(len(cForward)): 
            if len(cForward[f]) > 4: 
                fstrand1 = cForward[f][0].getPos1() 
                fstrand2 = cForward[f][0].getPos2() 
            
                for cReverse in reverse: 
                    for r in range(len(cReverse)): 
                        if (cReverse[r][0].getChromRef() != cForward[f][0].getChromRef()): break # break if chrom ref does not match
                        if len(cReverse[r]) > 4: 
                            rstrand1 = cReverse[r][0].getPos1() 
                            rstrand2 = cReverse[r][0].getPos2() 
                            if ( abs(fstrand1 - rstrand1) < 30000 and abs(fstrand2 - rstrand2) < 30000): 
                                fsize = len(cForward[f])
                                rsize = len(cReverse[r]) 
                                biggest = 0 
                                smallest = 0
                                # this is forthe regular breakpoints print 
                                for i in cForward[f]: 
                                    text = str(f) + "\t" + i.getID() +"\t"+str(i.getSamFlag()) +"\t"+ str(i.getChromRef()) + "\t" + str(i.getMapQual()) +"\t"+str(i.getPos1())+"\t" +str(i.getPos2())
                                    writeFile.write(text+"\n")
                                for i in cReverse[r]: 
                                    text = str(r) + "\t" + i.getID() +"\t"+str(i.getSamFlag()) +"\t"+ str(i.getChromRef()) + "\t" + str(i.getMapQual()) +"\t"+str(i.getPos1())+"\t" +str(i.getPos2())
                                    writeFile.write(text+"\n")                         
                                writeFile.write("\n")
                                break 

# prints out the summaries 
def printSum(forward, reverse, writeFile):
    for cForward in forward:
        for f in range(len(cForward)): 
            if len(cForward[f]) > 4: 
                fstrand1 = cForward[f][0].getPos1() 
                fstrand2 = cForward[f][0].getPos2() 
            
                for cReverse in reverse: 
                    for r in range(len(cReverse)): 
                        if (cReverse[r][0].getChromRef() != cForward[f][0].getChromRef()): break # break if chrom ref does not match
                        if len(cReverse[r]) > 4: 
                            rstrand1 = cReverse[r][0].getPos1() 
                            rstrand2 = cReverse[r][0].getPos2() 
                            if ( abs(fstrand1 - rstrand1) < 30000 and abs(fstrand2 - rstrand2) < 30000): 
                                fsize = len(cForward[f])
                                rsize = len(cReverse[r]) 
                                biggest = 0 
                                smallest = 0

                                # this is for the summary print
                                cForward[f].sort(key=lambda x: x.getPos1())
                                writeFile.write("{}\t{}\t{}\t{}\t".format(f, cForward[f][0].getSamFlag(), cForward[f][0].getChromRef(), cForward[f][0].getPos1()))
                                pos1 = cForward[f][-1].getPos1() 

                                cForward[f].sort(key=lambda x: x.getPos2())
                                pos2 = cForward[f][-1].getPos2() 
                                writeFile.write("{}\t{}\t{}\t{}\n".format(cForward[f][0].getPos2(), pos1, pos2, len(cForward[f])))

                                # printing the rev here
                                cReverse[r].sort(key=lambda x: x.getPos1())
                                pos1 = cReverse[r][-1].getPos1() 
                                writeFile.write("{}\t{}\t{}\t{}\t".format(r, cReverse[r][0].getSamFlag(), cReverse[r][0].getChromRef(), cReverse[r][0].getPos1()))

                                cReverse[r].sort(key=lambda x: x.getPos2())
                                pos2 = cReverse[r][-1].getPos2() 
                                writeFile.write("{}\t{}\t{}\t{}\n".format(cReverse[r][0].getPos2(), pos1, pos2, len(cReverse[r])))
                       
                                writeFile.write("\n")
                                break 

# helper function for merging, creates a dictionary that checks for twice occurence and count it 
# then it checks the number of twice against half of total to see if it should merge 
def duplicate1(clist1, clist2): 
    tmpdict = dict()
    tmpdict1 = dict() # this is to hold the class 
    tmplist = list()  
    count = 0 
    total = 0

    for i in clist1: 
        text = i.getID()
        # text = text.split(".")
        a = text
        tmpdict[a] = 1 
        tmpdict1[a] = i # getting the class here 

    for i in clist2: 
        text = i.getID()
        # text = text.split(".")
        a = text

        if not a in tmpdict: 
            tmpdict[a] = 1 
        else: 
            tmpdict[a] += 1 
        
        tmpdict1[a] = i

    for key, value in tmpdict.items(): # counting the number of 2s 
        if value == 2: 
            count += 1 
        total+= 1 

    if (total/2) < count: # check and submit the list 
        for key,value in tmpdict1.items():
            tmplist.append(value)
    else: 
        return tmplist # submits a list of 0 if total/2 is greater 
    return tmplist

# goes down the list of each chrom from the back to avoid mismatch. 
# returns a list of lists of classes 
def merging(superSet): 
    # merging into two different set 
    tmpList = list() 

    for flist in superSet: # this is getting a specific list for a chrom of the 4 lists inside
        for first in range(len(flist)-1, 0, -1): 
            for second in range(first-1, -1, -1): 
                randlist = duplicate1(flist[first], flist[second]) 
                if len(randlist) > 0: 
                    flist[second] = randlist
                    del flist[first] 
                    break
        tmpList.append(flist) 
    return tmpList

         
# the first half of main can be used as a function 
# this runs the entire program of the cluster algorithm 
def main(): 
    readFile = sys.argv[1] # samfile to be read 
    if (".sam" not in sys.argv[1]): 
        print("Not a SAM file/type")
        sys.exit()
    filename = readFile.find(".sam") 
    filenames = readFile[:filename]
    wFile = "good_" + filenames + "-result" # writing the name
    sFile = "summary_" + filenames + "-result" 

    writeFile = open(wFile, 'w') 
    sumFile = open(sFile, 'w') 
    for1 = [65, 129, 2161, 2225] #samflags to look out for
    rev1 = [113, 177, 2113, 2177] #samflags to look out for
    chrom = dict() 
    chrom_counter = 0
    forList = []
    revList = []

    with open(readFile, 'r') as f: 
        for line in f: 
            content = line.split('\t')

            if content[0] == "@SQ" : continue 
            if len(content) < 8: continue 
            if int(content[4]) < 20 or content[6] != "=": continue 

            entry = location(content[0], content[1], content[2], content[3],\
                            content[4], content[6], content[7])

            if abs(entry.getPos1() - entry.getPos2()) <= 1000000 or (entry.getPos1() - entry.getPos2()) < 0: 
                continue 
            
            # checks for one of the 4 chroms and give the counter object to it 
            if (entry.getChromRef() not in chrom): 
                chrom[entry.getChromRef()] = chrom_counter
                chrom_counter += 1
                forList.append(list())
                revList.append(list()) 

            counter = chrom[entry.getChromRef()] # get chrom position
            
            # checks to see if it is a forward or reverse reads and set the boolean for it 
            if entry.getSamFlag() in for1: 
                entry.convertFor() 
            elif entry.getSamFlag() in rev1: 
                entry.convertRev() 
            else: 
                continue 
            
            # checks to see if there are any clusters the entry matches and submit it there else 
            # return false and no list
            inCluster = False
            if (entry.getRev() == True): 
                """This uncommented section works in theory so we should be good to? """
                inCluster, revList[counter] = search(revList[counter], entry, 10000)
            elif(entry.getFor() == True): 
                inCluster, forList[counter] = search(forList[counter], entry, 10000) 
            # only if it is false, comes through here.
            # depending on forward or reverse strand send it to its respective chrom list 
            if inCluster == False: 
                tmpList = list() 
                tmpList.append(entry) 
                if(entry.getRev() == True): 
                    revList[counter].append(tmpList) 
                else: 
                    forList[counter].append(tmpList) 
        # merging duplicated clusters into one,
        # only exceptions are recurrence clusters that can safely avoid detection 
        forList = merging(forList)  
        revList = merging(revList) 
        printCluster(forList, revList, writeFile)
        printSum(forList, revList, sumFile)

    # making sure it writes by flushing it and closing the files. 
    writeFile.write("\n")
    writeFile.flush()
    writeFile.close()
    print("done with clustering") 
                          
if __name__ == "__main__": 
    main() 
