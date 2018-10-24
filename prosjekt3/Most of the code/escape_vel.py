from solarsystem import solsys
import numpy as np

ss = solsys()
###
g = 4*np.pi**2
v = np.sqrt(g)
pos = [1, 0, 0]
vel = [0, 0.99*np.sqrt(2)*v/365.2242, 0]
time = 2000
###
ss.addBody("EARTH", vel, pos, 1e-6)
ss.exportInitialValues("startValues.npy")
ss.simulate("startValues.npy", "orbits.txt", time = time, dt = 1e-4)