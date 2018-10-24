from solarsystem import solsys
import numpy as np

ss = solsys()
###
g = 4*np.pi**2
v = np.sqrt(g)
pos = [1, 0, 0]
vel = [0, 0.9*np.sqrt(2)*v/365.2242, 0]
time = 2000
###
ss.addBody("Testplanet", vel, pos, 1e-6)
ss.exportValues("c_test_values.npy")
ss.simulate("c_test_values.npy", "c_test_orbits.txt", time = time)