#!/usr/bin/python

import sys 
""" what : this is a post processing script for the outputs of reading3_1.cpp/cluster.py that searches for potential 
           breakpoints more accurately leaving only the recurrences and some other ones.

    input : python2.7 search_transposeable.py <name_of_SRR/RAL-name> (ex. ZI213)
    output : summary_<name>-breakpoints , summary_<name>-results 

    sidenote: it needs recurrence.csv, knownInversion.txt, notKnownInversion.txt, dmelMapTable.txt, TE_sites.txt

    #TODO: 3. Will keep both dmelmaptable and TE_sites.txt for further notice

""" 
# regular entry but mainly for TE 
class location:
    trans_size = 0 # not sure if list or class has a built in size check?  
    def __init__(self, chrom_ref, position1, position2):
        self.chrom_ref = chrom_ref
        self.position1 = position1 
        self.position2 = position2 
        location.trans_size += 1 

    def getposition1(self): 
        return self.position1 
    def getposition2(self): 
        return self.position2 
    def getchrom_ref(self):
        return self.chrom_ref 
    def getSize(self):
        return self.trans_size 
    def displayInfo(self): 
        print "this is the current size of this particlar class -chrom_ref " + self.chrom_ref + " -position1 " + str(self.position1) + " -position2 " + str(self.position2) + " " 

# this is for the breakpoint that comes in
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
        self.status = False
    
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
    def convertStatus(self): 
        self.status = True 
    def getStatus(self): 
        return self.status
    def displayInfo(self):
        print " this is the breakpoint here " + str(self.samflag) + " " + self.chrom_ref + " " + str(self.low_position1) + " " + str(self.low_position2) + " " + str(self.high_position1) + " " + str(self.high_position2) 

name = str(sys.argv).split(" ")  # would need to change this and allow multiple files to be read here

"""
# check out the tranposeable elements here and store them into a list or whatever that is easier to do so 
# then read in everything of the entire line - class - then compare only the first and second locations 
#                   ----> check to see if it is within the range tho.  
""" 
te_objects = []  # holds TE elements in class 
break_objects = []  # holds cluster in class 


# getting the value of the TE here 
with open("TE_sites.txt", "r") as f: 
    for line in f: 
        content = line.split() 
        decoy = location( content[0], int(content[3]), int(content[4]) )
        # print "this is the current counter after each line " + str(decoy.getSize())
        te_objects.append(decoy)   


# the idea is to store the breakpoints into their own separate list
# to make it easier to print and compare them.

break_store = [] # this will hold a list of classes 
check_sum = 0 
summaryResults = "summary_"+name[1].lstrip("'").rstrip("']")+"-result" 
with open(summaryResults, "r") as f: 
    for line in f: 
        # if there's an empty line as seen in the summary file, then it breaks into another list.
        if line in ['\n', '\r\n']:
            break_store.append(break_objects)
            break_objects = [] 
            check_sum += 1 
        else: 
            content = line.split() 
            decoy = b_location( int(content[0]), int(content[1]), content[2], int(content[3]), int(content[4]), int(content[5]), int(content[6]), int(content[7]))
            break_objects.append(decoy) 

# if there's only one cluster so it will never go through the other one
# because there's no newline 
if check_sum == 0 : 
    break_store.append(break_objects)

"""
# searches and compare the TE element with the cluster and see if it is close
# if one part matches then count, if at least like 50 - 80 % matches then add it to the list
"""
isTE = [] # clusters that did not pass the TE elements 
notTE = [] # clusters that did pass the TE elements 

counter = 0 # this is to keep track of the number of reads within TE at least 80% 
checklist = False 
for index, cur_list in enumerate(break_store): # getting each individual list

    in_transrange = False
    for t_index, cur_breakpoint in enumerate(cur_list): # iterating individually through chosen list

        in_transrange = False 
        for te in te_objects:
            if cur_breakpoint.getchrom_ref() == te.getchrom_ref(): # making sure the location is the same
            # somehow set_5_8 goes through. 
                if te.getposition1() <= cur_breakpoint.getlowposition1() <= te.getposition2() or te.getposition1() <= cur_breakpoint.gethighposition1() <= te.getposition2()\
                 or te.getposition1() <= cur_breakpoint.getlowposition2() <= te.getposition2() or te.getposition1() <= cur_breakpoint.gethighposition2() <= te.getposition2():
                    in_transrange = True
                    break
        if (in_transrange == True):
            isTE.append(cur_list)
            break 
    if in_transrange == False:
        notTE.append(cur_list)

# # TODO : what is this one tho? I can't recall?
# summaryBreakpoints = "summary_"+name[1].lstrip("'").rstrip("']")+"-breakpoints"
# with open(summaryBreakpoints, 'w') as f: 
#     for breakpoints in notTE: 
#         for index in breakpoints: 
#             f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format( index.getsamflag(), index.getchrom_ref(), index.getlowposition1(), index.getlowposition2(), index.gethighposition1(), index.gethighposition2(), index.getKnown() ))

# printing a new summary file here
with open(  ) # TODO HERE MAN 



goodReadFile = "good_"+name[1].lstrip("'").rstrip("']")+"-result" # ex. good_857-result
detailedBreakpoints = name[1].lstrip("'").rstrip("']")+"-breakpoints" # 857-result is the most important file in here
writeDetailedBreakpoints = open(detailedBreakpoints, 'w')

# this is opening good_<SRR>-results which is has the full info on the clusters
with open(goodReadFile) as f:
    for line in f:
        content = [x.strip() for x in line.split("\t")]
        if len(content) != 7: # to avoid getting an error for index out of range due to the empty line separating them? 
            continue 
        for i in newBreakStore: 
            for j in i: 
                if content[0] == str( j.getindex() ): 
                    writeDetailedBreakpoints.write(line)
                    break  
                    
                
