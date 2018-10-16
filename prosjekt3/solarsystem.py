import numpy as np

class solsys():

    def __init__(self):
        self.planets = []
        self.planPos = np.array()
        self.planVel = np.array()
        self.time = 0

    def values(self, valueFile):
        valfil = open(valieFile, "r")
        #do reading stuff
        #assign values to initiated variables

    def forward(self, intFunk, howFar, timeUnit): #runs simulation "howFar" "timeUnit"s forwards
        #something about converting time to seconds
        #case: seconds, hours, days, years
        #integrating "howFar" forwards using "intFunk"
            #sending start info to file
            #reading and running in julia
            #writing to file in julia
            #reading back
        #return/plot result