import os
"""This file will check for similar clusters throughout the files and
 other files to make sure there's no redundant everywhere
 
 Right now: run this later and do it after everything is working. 
 """ 
# to call this separately: 
# python similarclusters.py <# of SRR> 
# as long as there is 1+ clusters of almost identical reads, check for them all. 
class b_location: 
    def __init__(self, index, samflag, chrom_ref, low_position1, low_position2, high_position1, high_position2, size): 
        self.index = index 
        self.samflag = samflag 
        self.chrom_ref = chrom_ref 
        self.low_position1 = low_position1
        self.low_position2 = low_position2
        self.high_position1 = high_position1 
        self.high_position2 = high_position2 
        self.size = size 
    
    # class functions 
    def getsamflag(self): 
        return self.samflag
    def getlowposition1(self): 
        return self.low_position1
    def getlowposition2(self): 
        return self.low_position2
    def gethighposition1(self): 
        return self.high_position1
    def gethighposition2(self):
        return self.high_position2
    def getchrom_ref(self): 
        return self.chrom_ref 
    def getindex(self): 
        return self.index 
    def getSize(self):
        return self.size
    def displayInfo(self):
        print " this is the breakpoint here " + str(self.samflag) + " " + self.chrom_ref + " " + str(self.low_position1) + " " + str(self.low_position2) + " " + str(self.high_position1) + " " + str(self.high_position2) 

def main(): 
    arg_content = str(sys.argv).split() 
    content = arg_content[1].lstrip("['").rstrip(",']")
    lineNumber = int(content) # line number as in how many set files are there 
    fileposition = arg_content[2].lstrip("['").rstrip(",']") # fileposition = 
    original_breakpoints = []
    recurrence = []
    # reading every set file for the breakpoints
    for i in xrange(1,int(lineNumber)+1): 
        filename = "summary-set_{0}_{1}-breakpoiints".format(fileposition, i)
        readfile = open(filename, 'r')
        for i in readfile:
            content = [x.strip for x in i.split("\t")]
            decoy = b_location( int(content[0]), int(content[1]), content[2], int(content[3]), int(content[4]), int(content[5]), int(content[6]), int(content[7]))
            original_breakpoints.append(decoy) 
        readfile.close()
    
        


main()  