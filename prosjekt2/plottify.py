import matplotlib.pyplot as plt
import scipy as sc
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

n = []
counter = [] 
with open("rotated.txt") as f:
    content = f.readlines()
content = [(x.strip('\n').split()) for x in content]
for i in content:
    n.append(int(i[0]))
    counter.append(int(i[1]))

plt.figure()
plt.plot(n, counter)
plt.show()
plt.savefig("requiredRotations.pdf", bbox_inches = "tight")