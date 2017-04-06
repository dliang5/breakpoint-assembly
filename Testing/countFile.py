import subprocess as sp 
with open("SRRname", 'r') as f: 
    yes =  sum(1 for _ in f) 
    f.seek(0) 
    path = "./rename.sh " + str(1) + " " + str(yes)
    sp.Popen( path )
