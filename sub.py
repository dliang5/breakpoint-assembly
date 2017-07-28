with open("Breakpoints.csv", 'r') as read:
    bob = read.readline()
    for line in read:
        content = line.split(",")
        if content == "in between":
            continue
        print("the distance is {}".format(int(content[1]) - int(content[2])))
