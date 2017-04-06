import subprocess 
for i in xrange(0,20): 
    path = "echo " + str(i) 
    subprocess.Popen(path, shell=True, executable='/bin/bash') 
    subprocess.Popen("wait $!", shell=True, executable='/bin/bash')
    subprocess.Popen("echo bob " + str(i) , shell=True, executable='/bin/bash')
