# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 15:20:04 2018

@author: Bendik
"""

import seaborn
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, savefig, xscale, yscale

N = [10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
fx = [0.08811876090151421, 0.004634865424186327, 0.00043691327198678843, 4.3436014384114236e-5, 4.3410584003839165e-6, 4.3408078473925935e-7, 4.340766603464195e-8, 4.193607084543671e-9]


plot(N, fx)
xlabel("n"); ylabel("Relative error")
title("Plot of the relative error")
xscale("log"); yscale("log")
savefig("test_d.pdf")
show()