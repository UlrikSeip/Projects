# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 17:05:33 2018

@author: Bendik
"""

from numpy import zeros, exp, linspace, ones
from matplotlib.pyplot import plot, show, grid

n = 1000
a = ones(n-1)
a *= -1
b = ones(n)
b *= 2
c = ones(n-1)
c *= -1

B = zeros(n)
v = zeros(n)
b_ = zeros(n)
k = zeros(n)
B_ = zeros(n)
x = linspace(0, 1, n)

h = 1/(n+1)

def f(x):
    return (100*exp(-10*x))

u = 1 - (1 - exp(-10))*x - exp(-10*x)

for i in range(n):
    B[i] = h**2*f(x[i])

    
b_[0] = b[0]
B_[0] = B[0]
for i in range(1, n):
    temp = a[i-1]/b_[i-1]
    b_[i] = b[i] - temp*c[i-1]
    B_[i] = B[i] - temp*B_[i-1]

k[n-1] = B_[n-1]
i = n-2
while i >= 0:
    k[i] = B_[i] - c[i]*k[i+1]/b_[i]
    i-=1
    
print(k)
print(b_)

print(k/b_)

plot(x, k)
plot(x, u)
show()

