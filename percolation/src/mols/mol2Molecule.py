import re

def readLineAndSplit(line):
    return re.split("\s+",line.strip())

class Mol2Molecule:
    """
    Description of Mol2Molecule

    Attributes:
        triposAtom (type):
        triposBond (type):

    Args:
        mol2FileName (undefined):

    """
    
    def __init__(self,mol2FileName):
        self.triposAtom = []
        self.triposBond = []
        self.readMol2File(mol2FileName)
        
    def readMol2File(self,mol2FileName):        
        with open(mol2FileName) as mol2File:
            self.mol2File = mol2File
            for line in mol2File:
                self.readTriposMolecule()
                self.skipToNextTriposSection()
                self.readTriposAtom()
                self.skipToNextTriposSection()
                self.readTriposBond()
                
    def skipToNextTriposSection(self):
        line = self.mol2File.readline()
        while line[0] != "@":
            line = self.mol2File.readline()
            
    def readTriposMolecule(self):
        self.mol2File.readline()
        line = readLineAndSplit(self.mol2File.readline())
        self.numAtoms = int(line[0])
        self.numBonds = int(line[1])
        
    def readTriposAtom(self):
        for atom in range(self.numAtoms):
            line = readLineAndSplit(self.mol2File.readline())
            self.triposAtom.append((int(line[0]),float(line[2]),
            float(line[3]),float(line[4])))
            
    def readTriposBond(self):
        for bond in range(self.numBonds):
            line = readLineAndSplit(self.mol2File.readline())
            self.triposBond.append((int(line[1]),int(line[2])))
        
    