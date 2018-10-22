# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:07:42 2018

@author: Bendik
"""

import solarsystem


sSystem0 = solarsystem.solsys()
#sSystem0.addAllPlanets()
sSystem0.addBodyFromFile("EARTH")
sSystem0.addBodyFromFile("JUPITER")
#sSystem0.exportValues("noe.npy")
#sSystem0.simulate("testValues.npy", "orbitsTest.txt", time = 2)

sSystem1 = solarsystem.solsys()
sSystem1.addBodyFromFile("EARTH")
sSystem1.addBodyFromFile("JUPITER")
sSystem1.planets[1].mass *= 10


sSystem2 = solarsystem.solsys()
sSystem2.addBodyFromFile("EARTH")
sSystem2.addBodyFromFile("JUPITER")
sSystem2.planets[1].mass *= 1000

sSystem3 = solarsystem.solsys()
sSystem3.addBodyFromFile("EARTH")
sSystem3.addBodyFromFile("JUPITER")


sSystem4 = solarsystem.solsys()
sSystem4.addAllPlanets()

