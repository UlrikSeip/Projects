# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 14:27:17 2018

@author: Bendik
"""

import seaborn
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, savefig, xscale, yscale, legend

N = [10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
#b = [4.5039e-6, 6.6371e-6, 8.14618e-5, 0.00100219, 0.01926, 0.151202, 1.32222]
#c = [1.7383e-6, 9.9951e-6, 9.09433e-5, 0.00132338, 0.0171799, 0.148494, 1.13034]
#b = [2.12013e-6, 1.36296e-5, 0.000148096, 0.00137763, 0.019827, 0.15274, 1.32813]
#c = [1.75147e-6, 1.28264e-5, 0.000169745, 0.00110894, 0.0121225, 0.130161, 1.16427]
b = [1.44986e-6, 1.10973e-5, 9.99231e-5, 0.00144302, 0.0178639, 0.149568, 1.31139]
c = [1.12989e-6, 9.64741e-6, 0.000138801, 0.00159209, 0.0124083, 0.124903, 1.12878]
#with the average of 100 times

plot(N, b)
plot(N, c)
xlabel("n"); ylabel("Time (s)")
title("Plot of the time")
legend(["b", "c"])
xscale("log")#; yscale("log")
savefig("test_c.pdf")
show()
