#!bin/bash
""" This file will read in the Dmel line by line and parse the chrom region and 
the positions 1 and 2""" 

with open("DmelMapTable.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split(',')]
        """ this gets the chrom_region and positions together"""
        chrom_info = [x.strip() for x in content[2].split(":")]
        # print chrom_info
        """ this will now separate the positiosn from .. """
        if ( len(chrom_info) == 2):  
            pos_info = [x for x in chrom_info[1].split("..")]
            print chrom_info[0], pos_info[0], pos_info[1] 
        else: 
            continue 