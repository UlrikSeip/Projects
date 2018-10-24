import numpy as np
import solarsystem
import numpy as np

"""
sSystem0 = solarsystem.solsys()
sSystem0.addBodyFromFile("EARTH")
sSystem0.addBodyFromFile("JUPITER")
sSystem0.exportValues("kojr_e0.npy")
sSystem0.simulate("kojr_e0.npy", "ut_e0.txt", time = 30, dt = 1e-5, masses = np.array([sSystem0.planets[0].mass, sSystem0.planets[1].mass]), acc_func=2)


sSystem1 = solarsystem.solsys()
sSystem1.addBodyFromFile("EARTH")
sSystem1.addBodyFromFile("JUPITER")
sSystem1.planets[1].mass *= 10


sSystem2 = solarsystem.solsys()
sSystem2.addBodyFromFile("EARTH")
sSystem2.addBodyFromFile("JUPITER")
sSystem2.planets[1].mass *= 1000
"""

sSystem3 = solarsystem.solsys()
sSystem3.addBodyFromFile("EARTH")
sSystem3.addBodyFromFile("JUPITER")
sSystem3.addSun("SUN", 1)
sSystem3.exportValues("kojr_e0.npy")
sSystem3.simulate("kojr_e0.npy", "ut_e0.txt", time = 20, dt = 1e-5, masses = np.array([sSystem3.planets[0].mass, sSystem3.planets[1].mass, sSystem3.planets[2].mass]), acc_func=3)

"""
sSystem4 = solarsystem.solsys()
sSystem4.addAllPlanets()
sSystem4.addSun("Sun", 1)
sSystem4.exportValues("kojr_e0.npy")
sSystem4.simulate("kojr_e0.npy", "ut_e0.txt", time = 30, dt = 1e-5, masses = np.array([sSystem4.planets[0].mass, sSystem4.planets[1].mass, sSystem4.planets[2].mass]), acc_func=3)
"""

