import seaborn
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from celestialBodies import celestialBodies

class solsys():

    def __init__(self):
        self.time = 0
        self.planets = []

    def importValues(self, filename):
        with open(filename, "r") as infile:
            data = infile.readlines()
        counter = 0
        self.planPos = np.zeros((3, int(len(data)/3)))
        for i in range(len(data)):
            if counter == 3:
                counter = 0
            if counter == 0:
                self.planPos[0, int(np.ceil((i+1)/3)-1)] = float(data[i])
            if counter == 1:
                self.planPos[1, int(np.ceil((i+1)/3)-1)] = float(data[i])
            if counter == 2:
                self.planPos[2, int(np.ceil((i+1)/3)-1)] = float(data[i])
            counter += 1

    def addBody(self, name, vel0, pos0, mass):
        self.planets.append(celestialBodies(name, vel0, pos0, mass))
        
    def addSun(self, name, mass):
        pos0 = [0,0,0]
        vel0 = [0,0,0]
        for i in self.planets:
            vel0 += i.vel0*i.mass
        vel0 = vel0/mass
        self.planets.append(celestialBodies(name, vel0, pos0, mass))

    def plottXYOrbit(self):
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[0])
        plt.show()
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[1])
        plt.show()
        #plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[2])
        #plt.show()
        plt.plot(self.planPos[0], self.planPos[1])
        plt.show()
        #print(len(self.planPos[1]))

    def simulate(self, time, dt = 1e-5):
        subprocess.Popen(["julia", "integrator.jl"])

    def storeSystem(self, filename):
        nroPlanets = len(self.planets)
        vels = np.zeros((nroPlanets, 3))
        poss = np.zeros((nroPlanets, 3))
        for i in range(nroPlanets):
            vels[i] = self.planets[i].vel
            poss[i] = self.planets[i].pos

            
if __name__ == '__main__' :    
    solarsystem = solsys()
    solarsystem.importValues("orbits.txt")
    solarsystem.plottXYOrbit()
    #print(solarsystem.planPos[0, -2])
    