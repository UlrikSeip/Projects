# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 19:38:52 2018

@author: Bendik
"""

from subprocess import call


call(["julia.exe", "b_test.jl"])
#call(['python.exe', "solarsystem.py"])
exec(open("./solarsystem.py").read())
