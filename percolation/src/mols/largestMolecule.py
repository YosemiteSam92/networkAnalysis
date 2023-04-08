import numpy as np
import sys

class LargestMolecule:
    def __init__(self, molGraph):
        self.trajSnapshot = molGraph.trajSnapshot
        self.data = molGraph.bondRecord.data
        self.timeStep = molGraph.bondRecord.timeStep
        self.atomTypes = molGraph.atomTypes
        self.chemicalElements = self.data.chemicalElements
        self.numAtoms = molGraph.numAtoms
        self.initLargestMolecule(molGraph)
        self.initBoxCoordinates()
        
    def initLargestMolecule(self, molGraph):
        self.largestMolecule = molGraph.largestMolecule
        self.largestMolecule = np.sort(self.largestMolecule)
        self.largestMoleculeNumAtoms = np.count_nonzero(self.largestMolecule)
        
    def initBoxCoordinates(self):
        self.xlo = self.data.xlo
        self.xhi = self.data.xhi
        self.ylo = self.data.ylo
        self.yhi = self.data.yhi
        self.zlo = self.data.zlo
        self.zhi = self.data.zhi
        self.boxLen = self.data.boxLen
            
    def printToXyz(self, directoryForLargest, directoryForOthers):
        fileNameLargest = directoryForLargest + "/largest_" + str(self.timeStep) + ".xyz"
        fileNameOthers = directoryForOthers + "/others_" + str(self.timeStep) + ".xyz"
        with open(fileNameLargest, "w") as fileLargest:
            self.printPreamble(fileLargest, self.largestMoleculeNumAtoms)
            with open(fileNameOthers, "w") as fileOthers:
                self.printPreamble(
                    fileOthers, 
                    self.numAtoms - self.largestMoleculeNumAtoms
                )
                for atomId in range(1, self.numAtoms + 1):
                    self.printXyzLine(atomId, fileLargest, fileOthers)
    
    def printPreamble(self, fileHandle, numAtoms):
        print(numAtoms, file = fileHandle)
        print(
            self.xlo, self.xhi, 
            self.ylo, self.yhi, 
            self.zlo, self.zhi, 
            file = fileHandle
        )
                
    def printXyzLine(self, atomId, fileLargest, fileOthers):
        xyzLine = self.createXyzLine(atomId)
        if atomId in self.largestMolecule:
            print(xyzLine, file = fileLargest)
        else:
            print(xyzLine, file = fileOthers)
                        
    def createXyzLine(self, atomId):
        chemicalElement = self.getChemicalElement(atomId)
        atomX = self.getAtomXasString(atomId)
        return (
            chemicalElement + " " + atomX[0] + " " + atomX[1] + " " + atomX[2]
        )
            
    def getChemicalElement(self, atomId):
        return self.chemicalElements[self.atomTypes[atomId - 1] - 1]

    def getAtomXasString(self, atomId):
        scaledX = self.trajSnapshot.x[atomId]
        return (
            str(self.xlo + scaledX[0] * self.boxLen[0]),
            str(self.ylo + scaledX[1] * self.boxLen[1]),
            str(self.zlo + scaledX[2] * self.boxLen[2])
        )
