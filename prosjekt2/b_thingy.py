import numpy.linalg as npl
import numpy as np

n = 5
a = np.ones((n, n))

"""
for i in range(n):
    for j in range(n):
        if ((j == i+1) or j == i-1):
            a[i, j] = 2                     # b = 2
        elif (i == j):
            a[i, j] = 5                     # a = 5
"""
eigVal, eigVec = npl.eig(a)

print(eigVal)
print(eigVec)