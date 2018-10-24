import numpy as np
import solarsystem

"""
sSystem0 = solarsystem.solsys()
sSystem0.addBodyFromFile("EARTH")
sSystem0.addBodyFromFile("JUPITER")
sSystem0.exportValues("kojr_e0.npy")
sSystem0.simulate("kojr_e0.npy", "ut_e0.txt", time = 30, dt = 1e-3, acc_func=2)


sSystem1 = solarsystem.solsys()
sSystem1.addBodyFromFile("EARTH")
sSystem1.addBodyFromFile("JUPITER")
sSystem1.planets[1].mass = 10 * sSystem1.planets[1].mass
sSystem1.exportValues("kojr_e1.npy")
sSystem1.simulate("kojr_e1.npy", "ut_e1.txt", time = 30, dt = 1e-5, acc_func=2)


sSystem2 = solarsystem.solsys()
sSystem2.addBodyFromFile("EARTH")
sSystem2.addBodyFromFile("JUPITER")
sSystem2.planets[1].mass = 1000 * sSystem2.planets[1].mass 
sSystem2.exportValues("kojr_e2.npy")
sSystem2.simulate("kojr_e2.npy", "ut_e2.txt", time = 30, dt = 1e-5, acc_func=2)


sSystem3 = solarsystem.solsys()
sSystem3.addBodyFromFile("EARTH")
sSystem3.addBodyFromFile("JUPITER")
sSystem3.addSun("SUN", 1)
sSystem3.exportValues("kojr_e3.npy")
sSystem3.simulate("kojr_e3.npy", "ut_e3.txt", time = 20, dt = 1e-5, acc_func=3)

"""
sSystem4 = solarsystem.solsys()
sSystem4.addAllPlanets()
sSystem4.addSun("SUN", 1)
sSystem4.exportValues("kojr_e4.npy")
sSystem4.simulate("kojr_e4.npy", "ut_e4.txt", time = 500, dt = 1e-3, acc_func=3)


