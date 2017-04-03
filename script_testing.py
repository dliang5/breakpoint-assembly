with open("good_ZI213-result", 'r') as f: 
    for line in f: 
        if line in ['\n', '\r\n']:
            # this occurs if there is a newline or empty line 
            print 'no'
        else: 
            # this occurs if there is no empty line but contents to read. 
            content = line.split() 
            print content[3]