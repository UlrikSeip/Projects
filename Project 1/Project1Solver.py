import scipy as sc
import matplotlib.pyplot as plt
from numba import jit
from copy import copy

dig = sc.ones(4)*2
dig1 = sc.ones(4)

u = sc.zeros(4)
b = sc.zeros(4)

u = sc.zeros(len(b))
b = u.copy()
u = 1
print(b)

@jit
def compute(d, e):

    u = sc.zeros(len(d))
    b = u.copy()
    u = 1
    print(b)
        #for i in range 
