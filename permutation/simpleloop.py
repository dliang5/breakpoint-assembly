# from itertools import combinations as comb
 
# def permutationTest(a, b):
#     ab = a + b
#     Tobs = sum(a)
#     under = 0
#     for count, perm in enumerate(comb(ab, len(a)), 1):
#         print("{} {} {} vs {}".format(count, perm, Tobs, sum(perm)))
#         if sum(perm) <= Tobs:
#             under += 1
#     return under * 100. / count
 
# treatmentGroup = [85, 88, 75, 66, 25, 29, 83, 39, 97]
# controlGroup   = [68, 41, 10, 49, 16, 65, 32, 92, 28, 98]
# under = permutationTest(treatmentGroup, controlGroup)
# print("under=%.2f%%, over=%.2f%%" % (under, 100. - under))

# with open("dmel-gene.txt", 'r') as f: 
#     for line in f:
#         content = line.split( )
#         name = content[8].split(";")
#         if "CR" in name[1]:
#             print name[1]
with open("dmel-gene.txt", 'r') as gene:
    for line in gene:
        content = line.split()
        if(content[0] not in ['2R', '3R', '2L', 'X', '3L']): continue
        name =content[8].split(";")
        if "CR" in name[1] and "+" in content[6]: continue
        print(content[0], content[3], content[4], name[1], content[6])