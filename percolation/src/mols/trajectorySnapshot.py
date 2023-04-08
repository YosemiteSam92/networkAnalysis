#!/usr/bin/env python3

import numpy as np

def splitRecordEntry(l):
    line = l.split(" ")
    return int(line[0]), line[2], line[3], line[4]


# Position information of individual trajectory frame
class TrajectorySnapshot:
    def __init__(self, trajRecord):
        self.initAttributes(trajRecord)
        self.initX()
        self.populateX()

    def initAttributes(self, trajRecord):
        self.x = {}  # key = atomId; value = [x,y,z] of atomId
        self.xAsString = {} # key = atomId; value = ["x","y","z"] of atomId
        self.numAtoms = trajRecord.numEntries
        self.trajRecordEntries = trajRecord.entries

    def initX(self):
        for i in range(1, self.numAtoms + 1):
            self.x[i] = np.array([0.0, 0.0, 0.0])

    def populateX(self):
        for entry in self.trajRecordEntries:
            id, x, y, z = splitRecordEntry(entry)
            self.xAsString[id] = [x, y, z]
            self.x[id] = np.array([float(x), float(y), float(z)])
            
        
