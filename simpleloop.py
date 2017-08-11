# """random testing for geneSearch"""
# with open("Breakpoints.csv", 'r') as file:
#     header = file.readline()
#     print("this is rare inversion part")
#     for entry in file:
#         print(entry.strip("\n"))
#         if "In between" in entry:
#             break
#     print("\n\n\n\n\nthis is common inversion part")
#     for entry in file:
#         print(entry.strip("\n"))
#         if "Unique set" in entry:
#             break
class inversion:
    def __init__(self, chromRef,lpos1, lpos2, hpos1, hpos2, offset, highOffset):
        self.chromRef = chromRef
        self.lpos1 = int(lpos1)
        self.lpos2 = int(lpos2)
        self.hpos1 = int(hpos1)
        self.hpos2 = int(hpos2)
        self.offset = int(offset)
        self.highOffset = int(highOffset)

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
    def gethighoffset(self): 
        return (self.highOffset) 
    def getChromRef(self): 
        return (self.chromRef) 
    def display(self): 
        print ("chrom = {} lp1 = {} lp2 = {} hp1 = {} hp2 = {} of = {} hof = {}".format(self.chromRef, self.lpos1, self.lpos2, self.hpos1, self.hpos2, self.offset, self.highOffset))

inversionList = list()
# parsing the inversion file out
with open("Breakpoints.csv", 'r') as inversionEntry:
    # getting rare inversion breakpoints here
    for times in range(2):
        currentList = list() # getting the first list for rare inversion
        while True:
            entry1 = inversionEntry.readline()
            if entry1.split(",")[0] == "Inversion ": continue # getting rid of the Inversion line and restarting 
            if "In between" in entry1: break # moving on to the common inversions
            entry2 = inversionEntry.readline()
            if not entry2 or "in betweens" == entry1.split(",")[0] == "in betweens": break
            content1 = entry1.split(",")
            content2 = entry2.split(",")
            if content1[0] == "In between" or content1[0] == "Unique set": break 
            print("{}\n{}".format(entry1.strip("\n"), entry2.strip("\n")))
            chrom = (content1[0].split("("))[1].split(")")[0]
            print("chrom = {} lp = {} lp = {} hp = {} hp = {} of = {} hof = {}".format(chrom, content1[1], content1[2], content2[1], content2[2], content1[3], content2[3]))
            invClass = inversion(chrom, content1[1], content1[2], content2[1], content2[2], content1[3], content2[3] )
            currentList.append(invClass)
        inversionList.append(currentList)

for i in range(2):
    print("this is {}".format(i))
    for entry in inversionList[i]:
        entry.display()
