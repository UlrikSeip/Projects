import matplotlib.pyplot as plt
import scipy as sc
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

n = []
simTrans = []
time = []

with open("rotated.txt") as f:                      #Leser og deler opp filen
    content = f.readlines()
content = [(x.strip('\n').split()) for x in content]
for i in content:                                   #Stapper data fra filen over i plottbare arrays
    n.append(int(i[0]))
    simTrans.append(int(i[1]))
    time.append(float(i[2]))

plt.figure()                                        #Plottestuff for "n, simTrans"
plt.plot(n, sc.array(simTrans)/(sc.array(n)))
plt.plot(n, time)
plt.xlabel('Matrix dim')
plt.ylabel('Similarity transformations per n')
plt.title('Visulisation of required transformations for larger matrices')
plt.legend(('simTrans/n', 'computation time'))
plt.show()
plt.savefig("requiredRotations.pdf", bbox_inches = "tight")
