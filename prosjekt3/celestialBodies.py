import numpy as np

class celestialBodies():

    def __init__(self, name, vel0, pos0):
        """
        Takes initial arrays, and stores them.
        
        Arguments:
            name {string} -- [name of body]
            vel0 {ndarray} -- [initial velocity]
            pos0 {ndarray} -- [initial position]
        """
        self.name = name
        self.vel0 = vel0
        self.pos0 = pos0
        self.vel = vel0
        self.pos = pos0