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
class b_location: 
    def __init__(self, index, ID, samflag, chrom_ref, MQ, position1, position2 ): 
        self.index = index 
        self.ID = ID 
        self.samflag = samflag 
        self.chrom_ref = chrom_ref 
        self.MQ = MQ 
        self.position1 = position1 
        self.position2 = position2 
    
    # class functions 
    def getID(self): 
        return self.ID 
    def getsamflag(self): 
        return self.samflag
    def getposition1(self): 
        return self.position1
    def getposition2(self): 
        return self.position2
    def getchrom_ref(self): 
        return self.chrom_ref 
    def getMQ(self): 
        return self.MQ 
    def getindex(self): 
        return self.index 
    def displayInfo(self):
        print " this is the breakpoint here " + str(self.samflag) + " " + self.chrom_ref + " " + str(self.position1) + " " + str(self.position2)
    
# how this works : 
# python search_transposeable name[1] or the breakpointfile here
print 'Argument List:', str(sys.argv)
name = str(sys.argv).split(" ") 
write_name = name[1].lstrip("'").rstrip("']")+"-breakpoints" # this is getting the name of the file being written into with actual breakpoints that pass the transponseable element test 
break_point_name = "summary-"+name[1].lstrip("'").rstrip("']")+"-result" 
write_file = open(write_name, "w")  

# check out the tranposeable elements here and store them into a list or whatever that is easier to do so 
# then read in everything of the entire line - class - then compare only the first and second locations 
#                   ----> check to see if it is within the range tho.   
trans_objects = []
break_objects = []  
# getting the value of the TE here 
with open("TE_sites.txt", "r") as f: 
    for line in f: 
        content = line.split() 
        decoy = location( content[0], int(content[3]), int(content[4]) )
        # decoy.displayInfo()
        # print "this is the current counter after each line " + str(decoy.getSize())
        trans_objects.append(decoy)   
# print "testing one more time" 
# for i in trans_objects: 
#     i.displayInfo() 

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
            print "testing"
            print content[0] 
            print content[1]
            print content[2]
            print content[3]
            print content[4]
            print content[5]
            print content[6]
            decoy = b_location( int(content[0]), content[1], int(content[2]), content[3], int(content[4]), int(content[5]), int(content[6]))
            break_objects.append(decoy) 
            decoy.displayInfo() 

# if there's only one cluster so it will never go through the other one.
if check_sum == 0 : 
    break_store.append(break_objects)

# searches and compare the TE element with the cluster and see if it is close
# if one part matches then count, if at least like 50 - 80 % matches then add it to the list
nope_number = []
yes_number = []
print "bob" 
for rand_list in break_store: # getting each individual list
    print " Here " 
    f_count = 0 # count to check
    s_count = 0 # checking the second position for it as well 
    for t_index, entry in enumerate(rand_list): # iterating individually through chosen list
        entry.displayInfo() 
        for j in trans_objects: # iterating through TE now 
            if entry.getchrom_ref() == j.getchrom_ref(): # making sure the location is the same
                # gotta check if let's say position1 or positoin2 is in the TE range 
                if j.getposition1() <= entry.getposition1() <= j.getposition2(): 
                    # print j.getposition1() + " against " + entry.getposition1()
                    # print entry.getposition1() + " against " + j.getposition2() 
                    f_count += 1   
                elif j.getposition1() <= entry.getposition2() <= j.getposition2(): 
                    # print j.getposition1() + " against " + entry.getposition2()
                    # print entry.getposition2() + " against " + j.getposition2() 
                    s_count += 1  
                billy = len(rand_list)
                print "this is the length of rand_list here should be 3 " + str(billy)
                print str(t_index)  
                if (len(rand_list) * .8 ) <= f_count or (len(rand_list) * .8) <= s_count: # if count is greater than the size of the cluster, then it's TE. 
                    nope_number.append(rand_list) # done checking so we push
                    print "nope"
                    continue
        if t_index == len(rand_list)-1: # if not a TE, push to yes_number
            print "yes" 
            yes_number.append(rand_list) 
            
print "mango" 
for i in break_store: 
    for j in i: 
        j.displayInfo()
print "checking yes now" 
for i in yes_number: 
    for j in i: 
        j.displayInfo()
# writing the cluster that is not have any issues. 
for j in yes_number: 
    print "bob"
    for count, i in enumerate(j): 
        # bobby = str(i.getposition1)
        # joe = str(i.getposition2)
    # print str(count) + "\t" + str(i.getindex()) + "\t" + i.getID() + "\t" + str(i.getsamflag()) + "\t" + i.getchrom_ref() + "\t" + str(i.getposition1) + "\t" + str(i.getposition2) 
        write_file.write(str(i.getindex()) + "\t" + i.getID() + "\t" + str(i.getsamflag()) + "\t" + i.getchrom_ref() + "\t" + str(i.getposition1()) + "\t" + str(i.getposition2()) + "\n") 
    write_file.write("\n")
