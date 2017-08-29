# creating a distribution graph in python since I'm not too sure how the cluster system will work
import numpy as np
import scipy.stats as stats 
import pylab as pl

# with open("dist_file.txt", 'r') as f: 
#     for line in f:
#         content = line.split()
#         h = line.split()[1] = line.split()[1].lstrip("[").rstrip("]").split(",")
#         fit = stats.norm.pdf(h, np.mean(h), np.std(h))
#         pl.plot(h, fit, '-o')
#         pl.his(h, normed=True)
#         pl.show()
a = list()
b = list()
with open("dist_file.txt", 'r') as f:
    for line in f: 
        content = line.split()
        if content[0] == '0':
            a.append(int(content[1]))
        elif content[0] == '1':
            b.append(int(content[1]))
    a.sort()
    b.sort()
    hmean = np.mean(a)
    hstd= np.std(a)
    fit = stats.norm.pdf(a, hmean, hstd)
    pl.plot(a, fit, '-o')
    pl.hist(a, normed=True)
    pl.show()