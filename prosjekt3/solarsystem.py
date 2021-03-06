#import seaborn
import numpy as np
import os
import subprocess
import sys
from celestialBodies import celestialBodies

class solsys():

    def __init__(self):
        self.time = 0
        self.planets = []

    def importValues(self, filename):
        with open(filename, "r") as infile:
            data = infile.readlines()
        counter = 0
        nroPlanets = len(self.planets)
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
    
    #not sure if this should be in the constructor or not
    #should be called after all the planets
    def addSun(self, name, mass):
        pos0 = [0,0,0]
        vel0 = [0,0,0]
        for i in self.planets:
#            pos0 += i.pos0*i.mass
            vel0 += i.vel0*i.mass
        #pos0 = pos0/mass
        vel0 = -vel0/mass
        print(pos0)
        self.planets.append(celestialBodies(name, vel0, pos0, mass))

    def plottXYOrbit(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D as plt3
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[0])
        plt.show()
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[1])
        plt.show()
        #plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[2])
        #plt.show()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for i in range(len(self.planets)):
            ax.plot(self.planPos[0], self.planPos[1], self.planPos[2])
        plt.show()
        #print(len(self.planPos[1]))

    def simulate(self, inFile, outFile, time = 10, dt = 1e-5, plott = "true", intgrat = 1, acc_func=1):
        names = []
        masses = []
        for i in range(len(self.planets)):
            names.append(self.planets[i].name)
            masses.append(self.planets[i].mass)            
        if intgrat == 1:
            simulation = subprocess.Popen(["julia", "velocity_verlet.jl", inFile, outFile, str(time), str(dt), plott, str(masses), str(names), str(acc_func)])
        elif intgrat == 2:
            simulation = subprocess.Popen(["julia", "forward_euler.jl", inFile, outFile, str(time), str(dt), plott, str(masses), str(names), str(acc_func)])
        simulation.wait() #waits for simulation to finish before doing anything else
        
    def addBodyFromFile(self, name):
        """
        Finds the planet 'name' in lesSolarsysteminfo.txt, and sets vel0, pos0 and mass to the subsequent lines in the file.
        'name' needs to be capitalised
        """
        with open("lesSolarsysteminfo.txt", "r") as infile:
            data = infile.readlines()
        k = 0
        for i in data:
            i = i.strip()
            if i == name:
                k = 1
            elif k == 1:
                mass = float(i)/2e30
                k = 2
            elif k == 2:
                pos0 = i.split(',')
                pos0[0] = float(pos0[0])
                pos0[1] = float(pos0[1])
                pos0[2] = float(pos0[2])
                pos0 = np.asarray(pos0)
                k = 3
            elif k == 3:
                vel0 = i.split(',')
                vel0[0] = float(vel0[0])
                vel0[1] = float(vel0[1])
                vel0[2] = float(vel0[2])
                vel0 = np.asarray(vel0)
                break
            
        if k == 0:
            print("Couldn't find " + name + " in the file. Make sure you wrote the name correctly and that it's capitalised")
            sys.exit()
        else:
            self.planets.append(celestialBodies(name, vel0, pos0, mass))

    def addAllPlanets(self):
        self.addBodyFromFile("MERCURY")
        self.addBodyFromFile("VENUS")
        self.addBodyFromFile("EARTH")
        self.addBodyFromFile("MARS")
        self.addBodyFromFile("JUPITER")
        self.addBodyFromFile("SATURN")
        self.addBodyFromFile("URANUS")
        self.addBodyFromFile("NEPTUNE")
        self.addBodyFromFile("PLUTO")

    def exportValues(self, filename): #.npy
        nroPlanets = len(self.planets)
        print(self.planets)
        vels = np.zeros(nroPlanets * 3)
        poss = np.zeros(nroPlanets * 3)
        print(nroPlanets)
        for i in range(0, nroPlanets*3, 3):
            vels[i] = self.planets[int(i/3)].vel[0]
            vels[i+1] = self.planets[int(i/3)].vel[1]
            vels[i+2] = self.planets[int(i/3)].vel[2]
            poss[i] = self.planets[int(i/3)].pos[0]
            poss[i+1] = self.planets[int(i/3)].pos[1]
            poss[i+2] = self.planets[int(i/3)].pos[2]
            #poss[i] = self.planets[i].pos
        """    
        with open(filename, 'w') as the_file:
            for i in range(nroPlanets):
                the_file.write(str(self.planets[i].pos)+'\n')
                the_file.write(str(self.planets[i].vel)+'\n')
        """
        #print(poss)
        np.save(filename, [poss, vels])


if __name__ == '__main__' :    
    solarsystem = solsys()
    #solarsystem.addAllPlanets()
    #solarsystem.addBodyFromFile("MERCURY")
    #solarsystem.addBodyFromFile("VENUS")
    #solarsystem.addBodyFromFile("EARTH")
    solarsystem.exportValues("testValues.npy")
    #solarsystem.simulate("testValues.npy", "orbitsTest.txt", time = 100) #time in days, can allso take dt, masses and names
    #solarsystem.importValues("orbitsTest.txt")
    #solarsystem.plottXYOrbit()
    #print(solarsystem.planPos[0, -2])




