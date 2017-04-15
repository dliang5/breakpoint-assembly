
with open("testingTE.txt", 'r') as f: 
    for i in f: 
        content = [x.strip() for x in i.split("\t")]
        for j in content: 
            print j