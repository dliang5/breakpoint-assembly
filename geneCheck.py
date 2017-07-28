#!/usr/bin/env python3
# David Liang dliang5 
# this file will check if the breakpoints interrupt a gene
# if not, it will find the closest gene to it

""" There will be two options of running this program 
    1. python geneCheck.py < [breakpoint] 
    2. python geneCheck.py [breakpoint] 

    Questions/Assumptions: 
    1. I will be using the two contigs and their length to match up
    2. First, start off by finding the start codon, if there's a stop codon check and match against it.
        2.a. Once a stop or start is found, gather all of the bases and create the sequence. If the string does not match 
        the contigs and what not, then move on.
    3. To increase accurate and the time complexities, only check for certain chromosomes of the contigs and the section of the breakpoints. 
    4. to make it automatic, it can take in a list of chrom and a list of breakpoints
        4.a. it reads them corresponding to each other
"""
""" 
In case I gotta look at the reverse strand as well
"""
def reverseComplement(sequence): 
    reverseSequence = "" 
    for nuc in sequence: 
        nucList = {"A":"T", "C":"G", "G":"C", "T":"A"}
        reverseSequence = nucList[nuc] + reverseSequence 
    return (reverseSequence)
        

def main(): 
    import sys, subprocess as sp, argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-ch', '--chrom', nargs='+', help='set chromosomes to check', required=True)
    parser.add_argument('-b', '--break', nargs='+', help='set breakpoints to read', required=True) 
    parser.add_argument('f', '--fasta', nargs="+", help="corresponding SAM file to read", required=True) 
    args = vars(parser.parse_args())

    contigTitle = list() 
    contigDict = dict() 
    # making the assumption that there is only one chrom and breakpoint to read to make it easier 
    # assumption, the details of the region is location with the name of the contigs
    with open(args["break"][0], 'r') as f:
        currentTitle = "" 
        for entry in f: 
            content = entry.split(" ") 
            if content[0] == ">":
                # what I'm getting here: name of contig, chrom, every positions, and finally length of the contigs
                tempTitle.append((content[0],args["chrom"][0], content[1], content[2], content[3], content[4]))
                currentTitle = content[0] 

            # getting the length of the contigs here until ">" 
            if len(content) == 1: 
                contigDict[currentTile] = contigDict[currentTitle] + len(content[0])
        
    # moving on... to reading the SAM or is it fasta file and producing the orfs
    # reading a file twice
    with open(args["break"]+"-gene", 'w') as f: 
        for times in range(2):
            myReader = sequenceAnalysis.FastAreader(args["fasta"])
            for head, seq in myReader.readFasta(): 
                
    
if __name__ == "__main__": 
    main() 