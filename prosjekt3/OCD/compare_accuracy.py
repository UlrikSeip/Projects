from solarsystem import solsys
import numpy as np

solarsystem = solsys()
solarsystem.addBodyFromFile("EARTH")
solarsystem.exportInitialValues("startValues.npy")
solarsystem.simulate("startValues.npy", "orbitsTest.txt", time = 10000, acc_func=1, int_func=1, dt = 1e-2)
solarsystem.simulate("startValues.npy", "orbitsTest.txt", time = 10000, acc_func=1, int_func=2, dt = 1e-2)