"""random testing for geneSearch"""
with open("Breakpoints.csv", 'r') as file:
    header = file.readline()
    print("this is rare inversion part")
    for entry in file:
        print(entry.strip("\n"))
        if "In between" in entry:
            break
    print("\n\n\n\n\nthis is common inversion part")
    for entry in file:
        print(entry.strip("\n"))
        if "Unique set" in entry:
            break