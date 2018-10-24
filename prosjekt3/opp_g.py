# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:56:29 2018

@author: Bendik
"""

import numpy as np
import solarsystem



sSystem3 = solarsystem.solsys()
sSystem3.addBodyFromFile("MERCURY")
sSystem3.planets[0].pos0 = [0.3075,0,0]
sSystem3.planets[0].vel0 = [0,-12.44,0]
sSystem3.exportValues("kojr_g.npy")
sSystem3.simulate("kojr_g.npy", "ut_e3.txt", time = 100, dt = 1e-3, acc_func=5)


