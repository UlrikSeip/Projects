import seaborn
import numpy as np
import matplotlib.pyplot as plt

class solsys():

    def __init__(self):
        self.time = 0

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

    def plottXYOrbit(self):
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[0])
        plt.show()
        plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[1])
        plt.show()
        #plt.plot(np.linspace(0, 1, len(self.planPos[1])),self.planPos[2])
        #plt.show()
        plt.plot(self.planPos[0], self.planPos[1])
        plt.show()
        print(len(self.planPos[1]))

solarsystem = solsys()
solarsystem.importValues("orbits.txt")
solarsystem.plottXYOrbit()
#print(solarsystem.planPos[0, -2])