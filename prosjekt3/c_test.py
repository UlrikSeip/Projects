from solarsystem import solsys
import numpy as np

ss = solsys()
g = 4*np.pi**2
v = np.sqrt(g)
pos = [1, 0, 0]
vel = [0, v/365.2242, 0]
time = 10
dt_values = [1e-1, 1e-2, 1e-3, 1e-4]
ss.addBody("Testplanet", vel, pos, 1e-6)
ss.exportValues("c_test_values.npy")
for i in dt_values: 
    ss.simulate("c_test_values.npy", "c_test_orbits.txt", time = time, dt = i)
    ss.importValues("c_test_orbits.txt")
    ss.plottXYOrbit()
