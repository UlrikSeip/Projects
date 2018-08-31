# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 17:05:33 2018

@author: Bendik
"""

from numpy import zeros, exp, linspace
import matplotlib.pyplot as plt

def f(x):
    return (100*exp(-10*x))

def calc(x, h):
    for i in range(n-1):
        a[i] = -1; c[i] = -1
    for i in range(n):
        b[i] = 2
        B[i] = h**2*f(x[i])

        
    b_[0] = b[0]
    B_[0] = B[0]
    for i in range(1, n):
        b_[i] = b[i] - a[i-1]*c[i-1]/b_[i-1]
        B_[i] = B[i] - a[i-1]*B_[i-1]/b_[i-1]

    k[n-1] = B_[n-1]
    i = n-2
    while i >= 0:
        k[i] = B_[i] - c[i]*k[i+1]/b_[i]
        i-=1
    return()

n = int(1e7)
a = zeros(n-1)
b = zeros(n)
c = zeros(n-1)

B = zeros(n)
v = zeros(n)
b_ = zeros(n)
k = zeros(n)
B_ = zeros(n)
x = linspace(0, 1, n)

h = 1/(n+1)

print(k)
print(b_)

print(k/b_)

plt.plot(x, k)
plt.show()
