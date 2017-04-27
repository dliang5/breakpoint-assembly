# this is a script for the post processing phase of reading.cpp or the file that searches for potential 
# inversion breakpoints. 

# What it should do, I guess it opens up the transposeable element from the reference genome 
# then it takes in the list of the potential inversion breakpoints and compare them. 
# -----> it any of them matches, throw that inversion breakpoint out and keep everything else into another file with actual 
#       -> points here. 

# how to use : python search_transposeable.py <name_of_file> (ex. ZI213) or $NAME  
# with script : qsub checking.sh ZI213 -> holds the variable of ZI213 and checks for that instead. 
#!/usr/bin/python

import sys 

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

# more specific for clusters over TE. 
# here's the template of the SRR reads 
# 194 65 2L 13604007 23063123 136006061 23063088
# basically, this shows the <cluster_index> <samflag> <chrom_region> <low_position1> <low_position2> <high_position1> <high_position2> 
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
"""
# how this works : 
# python search_transposeable name[1] or the breakpointfile here
it calls reads summary-ZI213-result 
"""
print 'Argument List:', str(sys.argv)
name = str(sys.argv).split(" ") 
summaryBreakpoints = "summary-"+name[1].lstrip("'").rstrip("']")+"-breakpoints" # this is getting the name of the file being written into with actual breakpoints that pass the transponseable element test 
break_point_name = "summary-"+name[1].lstrip("'").rstrip("']")+"-result" 
writeSummary = open(summaryBreakpoints, "w")  
"""
# check out the tranposeable elements here and store them into a list or whatever that is easier to do so 
# then read in everything of the entire line - class - then compare only the first and second locations 
#                   ----> check to see if it is within the range tho.  
""" 
trans_objects = []
break_objects = []  
# getting the value of the TE here 
with open("TE_sites.txt", "r") as f: 
    for line in f: 
        content = line.split() 
        decoy = location( content[0], int(content[3]), int(content[4]) )
        # print "this is the current counter after each line " + str(decoy.getSize())
        trans_objects.append(decoy)   


# getting the values of the clusters here from summary_*.txt  
break_store = [] # this will hold a list of list
# b_location ( index chrom_ref, posiiton1, position2, )
check_sum = 0 
with open(break_point_name, "r") as f: 
    for line in f: 
        if line in ['\n', '\r\n']:
            break_store.append(break_objects)
            #del break_objects[:] # deleting everything in break_objects. 
            break_objects = [] # recreating another list or emptying the current one
            check_sum += 1 
        else: 
            content = line.split() 
            decoy = b_location( int(content[0]), int(content[1]), content[2], int(content[3]), int(content[4]), int(content[5]), int(content[6]), int(content[7]))
            break_objects.append(decoy) 
            # decoy.displayInfo() 

# if there's only one cluster so it will never go through the other one.
if check_sum == 0 : 
    break_store.append(break_objects)

"""
reading the Dmel table to set up for the comparsion 
"""
dmel_holder = []
with open("DmelMapTable.txt", 'r') as f: 
    for i in f: 
        temp_dmel_store = []
        content = [x.strip() for x in i.split(',')]
        """ this gets the chrom_region and positions together"""
        chrom_info = [x.strip() for x in content[2].split(":")]
        # print chrom_info
        """ this will now separate the positiosn from .. """
        if ( len(chrom_info) == 2):  
            pos_info = [x for x in chrom_info[1].split("..")]
            temp_dmel_store.append( chrom_info[0] ) 
            temp_dmel_store.append( (pos_info[0]) ) 
            temp_dmel_store.append( pos_info[1] ) 
            dmel_holder.append(temp_dmel_store) 
            # print temp_dmel_store
        else: 
            continue 
# for i in dmel_holder: 
#     print i[0], i[1], i[2]  
"""check the values in yes_number to the dmel table to further more validate the results """ 
"""-----------------------------------------------------------------------------------==="""
"""
# searches and compare the TE element with the cluster and see if it is close
# if one part matches then count, if at least like 50 - 80 % matches then add it to the list
"""
transMatched = []
goodBreakpoints = []
sumBreakpoints = [] # this is me being lazy and appending the clusters that was checked by dmel table. 
print "bob" 
check_counter = 0 # in case if one matches a transposeable element, remove the corresponding cluster as it's kind of useless now. 
counter = 0 # this is to keep track of the number of reads within TE at least 80% 
checklist = False 
for index, rand_list in enumerate(break_store): # getting each individual list
    for t_index, entry in enumerate(rand_list): # iterating individually through chosen list
        in_transrange = False 
        for trans in trans_objects: # iterating through TE now 
            if entry.getchrom_ref() == trans.getchrom_ref(): # making sure the location is the same
            # somehow set_5_8 goes through. 
                if trans.getposition1() <= entry.getlowposition1() <= trans.getposition2() or trans.getposition1() <= entry.gethighposition1() <= trans.getposition2()\
                 or trans.getposition1() <= entry.getlowposition2() <= trans.getposition2() or trans.getposition1() <= entry.gethighposition2() <= trans.getposition2():
                    counter+=1
                if(counter >= entry.getSize()*.8 ):
                    cout << entry.displayInfo() << endl ; 
                    in_transrange = True 
        if (in_transrange == True):
            transMatched.append(rand_list)
            break 
    if in_transrange == False:
        goodBreakpoints.append(rand_list)


# compare dmel table to potential breakpoints
temp_dmel_holder = dmel_holder[:]
for t_index,j in enumerate(goodBreakpoints): 
    for count, i in enumerate(j):
        # print i.displayInfo()
        for index, dmel_line in enumerate(dmel_holder): 
            if (str(dmel_line[0]) == i.getchrom_ref()):
                if (int(dmel_line[1]) <= i.getlowposition1() <= int(dmel_line[2]) or int(dmel_line[1]) <= i.gethighposition1() <= int(dmel_line[2])):
                    # then check for the other end to see if they are around the same 
                    # temp_dmel_holder = dmel_holder[index:]
                    for temp_line in temp_dmel_holder:
                        if (str(temp_line[0]) == i.getchrom_ref()):
                            if (int(temp_line[1]) <= i.getlowposition2() <= int(temp_line[2]) or int(temp_line[1]) <= i.gethighposition2() <= int(temp_line[2])):
                            # if (i.getlowposition1() <= int(dmel_line[1]) <= i.gethighposition1() or i.getlowposition2() <= int(dmel_line[1]) <= i.gethighposition2()):
                                writeSummary.write(str(i.getindex()) + "\t" + str(i.getsamflag()) + "\t" + i.getchrom_ref() + "\t" + str(i.getlowposition1()) + "\t" + str(i.getlowposition2()) + "\t" + str(i.gethighposition1()) + "\t" + str(i.gethighposition2())+ "\t" + str(i.getSize()) + '\n' ) 
                                break # done with this position if everything matches 

writeSummary.close() 
               
recurrenceObjects = [] 
# getting the recurrence reads here . 
with open("recurrence.csv", 'r') as recurFile: 
    for entry in recurFile:
        content = [x.strip() for x in entry.split(" ,")]
        decoy = location(content[0], int(content[1]), int(content[2]))
        recurrenceObjects.append(decoy) 

# comparing the recurrence to the current potential breakpoints before going any further. 
with open(summaryBreakpoints, 'r') as f:  
    for i in f: 
        temp_sumBreakpoints = [] 
        content = [x.strip() for x in i.split("\t")]
        decoy = b_location( int(content[0]), int(content[1]), content[2], int(content[3]), int(content[4]), int(content[5]), int(content[6]), int(content[7]))
        for recurEntry in recurrenceObjects: 
            if decoy.getchrom_ref() != recurEntry.getchrom_ref(): 
                continue
            else: 
                if( (recurEntry.getposition1-1000 <= decoy.getlowposition1() <= recurEntry.getposition2+1000) or (recurEntry.getposition1-1000 <= decoy.getlowposition2() <= recurEntry.getposition2+1000 ):
                    if( (recurEntry.getposition1-1000 <= decoy.gethighposition1() <= recurEntry.getposition2+1000) or (recurEntry.getposition1-1000 <= decoy.gethighposition2() <= recurEntry.getposition2+1000 ):
                        continue 
                    else:       
                        temp_sumBreakpoints.append(decoy) 
                        sumBreakpoints.append(temp_sumBreakpoints)
        

""" writing a file of the actual reads based on the exisiting summary clusters
this will also check for the Dmel or the theoretical breakpoints here as well 
"""
read_file = "good_"+name[1].lstrip("'").rstrip("']")+"-result"
new_write = name[1].lstrip("'").rstrip("']")+"-breakpoints"
writing = open(new_write, 'w')

# this is opening good_<SRR>-results which is has the full info on the clusters
# why is it printing 10 of each element in a cluster? 
with open(read_file) as f:
    for line in f:
        content = [x.strip() for x in line.split("\t")]
        if len(content) != 7: # to avoid getting an error for index out of range due to the empty line separating them? 
            continue 
        #content = line.split("\t")
        for i in sumBreakpoints: 
            for j in i: 
                if content[0] == str( j.getindex() ): 
                    # start writing it 
                    writing.write(line)
		    writing.write("\n") 
                    break  
                    
                
