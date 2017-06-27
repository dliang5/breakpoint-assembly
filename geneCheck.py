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


def main(): 
    import sys, subprocess as sp, argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-ch', '--chrom', nargs='+', help='set chromosomes to check', required=True)
    parser.add_argument('-b', '--break', nargs='+', help='set breakpoints to read', required=True) 
    args = vars(parser.parse_args())

    # making the assumption that there is only one chrom and breakpoint to read to make it easier 
    
if __name__ == "__main__": 
    main() 