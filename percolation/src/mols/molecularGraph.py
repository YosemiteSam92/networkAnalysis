import numpy as np
import sys
from graph.graph_structs import MetaGraph, EdgeTag, EdgeEntry


def splitRecordEntry(l):
    line = l.split(" ")
    return int(line[1]), int(line[2]), int(line[3])


class MolecularGraph:
    def __init__(self, bondRecord, trajSnapshot):
        self.bondRecord = bondRecord
        self.trajSnapshot = trajSnapshot  # if no traj needed, just pass []
        self.numBondEntries = bondRecord.numEntries
        self.percolatingClusterFlag = False
        self.initAttributesFromData()
        self.initGraphAttributes()
        self.buildMolecularGraph()

    def initAttributesFromData(self):
        self.numAtoms = self.bondRecord.data.numAtoms
        self.atomicMasses = self.bondRecord.data.atomicMasses
        self.atomTypes = self.bondRecord.data.atomTypes
        self.totalMass = self.bondRecord.data.totalMass

    def initGraphAttributes(self):
        self.initBondAtomAndType()

    def initBondAtomAndType(self):
        self.bond_atom = {}  # key: atom id-1, value: list of atoms bonded to id
        self.bond_type = {}  # key: atom id-1, value: bond types of bonds involving id
        for i in range(1, self.numAtoms + 1):
            self.bond_atom[i] = []
            self.bond_type[i] = []

    def buildMolecularGraph(self):
        # vertex = atom; neighbors of vertex = bonding neighbors of atom
        for i in range(self.numBondEntries):
            btype, id1, id2 = splitRecordEntry(self.bondRecord.entries[i])
            self.createBond(btype, id1, id2)

    def createBond(self, btype, id1, id2):
        self.bond_atom[id1].append(id2)
        self.bond_atom[id2].append(id1)
        self.bond_type[id1].append(btype)
        self.bond_type[id2].append(btype)

    # ==== analyze molecular graph ==== #

    # COMPUTE MOLECULAR MASSES:
    # graph theory approach: count number of connected components of the
    # molecular graph using a depth first search (DFS), computing mass of
    # each component
    #
    # COMPUTE NUMBER OF INTRAMOLECULAR CROSSLINKS:
    # this is equal to the num of bonds connecting a C or N atom to a
    # previously seen C or N atom

    def computeMolecularMasses_numOfIntramolCrosslinks_largestMolecule(self):
        self.initDFSattributes()
        for atomId in range(1, self.numAtoms + 1):
            self.exploreMoleculeStartingFromUnseenAtom(atomId)
        self.handleLastMolecule()

    def initDFSattributes(self):
        self.visitedAtoms = np.zeros(self.numAtoms, dtype=int)
        self.numMolecules = 0
        self.numIntramolCrosslinks = 0
        self.molecularMasses = []
        self.currentMolecule = np.zeros(self.numAtoms, dtype=int)
        self.largestMoleculeMass = 0.0

    def exploreMoleculeStartingFromUnseenAtom(self, atomId):
        if self.visitedAtoms[atomId - 1]:
            return
        # new molecule found
        self.molecularMasses.append(0.0)
        self.numMolecules += 1
        self.updateLargestMolecule()
        self.currentMoleculeAtomIndex = 0
        self.depthFirstSearch(atomId, atomId)  # entry atom in a new molecule is
        # parent of itself

    def updateLargestMolecule(self):
        if self.numMolecules > 1:
            currentMoleculeMass = self.molecularMasses[self.numMolecules - 2]
            if currentMoleculeMass > self.largestMoleculeMass:
                self.updateLargestMoleculeAtoms()

    def updateLargestMoleculeAtoms(self):
        self.largestMolecule = self.currentMolecule
        self.largestMoleculeMass = self.molecularMasses[self.numMolecules - 2]
        self.currentMolecule = np.zeros(self.numAtoms, dtype=int)

    def handleLastMolecule(self):
        self.numMolecules += 1
        self.updateLargestMolecule()
        self.numMolecules -= 1

    def depthFirstSearch(self, startAtom, parentAtom):
        self.visitedAtoms[startAtom - 1] = 1
        self.addAtomToCurrentMolecule(startAtom)
        self.addAtomicMassToMolecularMass(startAtom)
        myNeighbors = self.bond_atom[startAtom]
        for bondAtom in myNeighbors:
            self.exploreUnseenNeighbor(bondAtom, startAtom, parentAtom)

    def addAtomToCurrentMolecule(self, atomId):
        self.currentMolecule[self.currentMoleculeAtomIndex] = atomId
        self.currentMoleculeAtomIndex += 1

    def addAtomicMassToMolecularMass(self, startAtom):
        atomicMass = self.atomicMasses[self.atomTypes[startAtom - 1] - 1]
        self.molecularMasses[self.numMolecules - 1] += atomicMass

    def exploreUnseenNeighbor(self, bondAtom, startAtom, parentAtom):
        if bondAtom == parentAtom:  # do not go back
            return
        if self.visitedAtoms[bondAtom - 1]:
            if self.thisNeighborFormsAnIntramolCrosslink(bondAtom, startAtom):
                self.numIntramolCrosslinks += 1
            return
        self.depthFirstSearch(bondAtom, startAtom)

    def thisNeighborFormsAnIntramolCrosslink(self, bondAtom, startAtom):
        # either bond atom or start atom is a nitrogen atom
        return self.atomTypes[bondAtom - 1] == 8 or self.atomTypes[startAtom - 1] == 8

    # PERCOLATING MOLECULE
    # find out whether a molecule closes a loop with itself, crossing the pbcs
    # in all dimensions - this algorithm suggested by Kevin Hoellring.
    # A pbc component containing atom i is composed of all those atoms that can
    # be reached from i via bonds that do not cross any pb.

    def getAtomX(self, atomId):
        return self.trajSnapshot.x[atomId]

    def isAnyMoleculePercolating(self):
        self.findAllPbcComponents()
        self.buildGraphAmongPbcComponents()
        self.findPercolatingMolecules()

    def findAllPbcComponents(self):
        self.initAttributesPbcComponents()
        for startAtom in range(1, self.numAtoms + 1):
            self.findPbcComponentContainingStartAtom(startAtom)

    def initAttributesPbcComponents(self):
        self.pbcDistanceThreshold = 0.9 # remember: dump file coords
                                        # are between 0 and 1
        self.visitedAtoms = np.zeros(self.numAtoms, dtype=int)
        self.numOfPbcComponents = 0
        self.atomsInPbcComponent = {}
        self.pbcComponentOfAtom = np.zeros(self.numAtoms, dtype=int)

    def findPbcComponentContainingStartAtom(self, startAtom):
        if not self.visitedAtoms[startAtom - 1]:
            self.addNewPbcComponent()
            self.depthFirstSearchPbcComponents(startAtom, startAtom)

    def addNewPbcComponent(self):
        self.numOfPbcComponents += 1
        self.atomsInPbcComponent[self.numOfPbcComponents - 1] = []

    def depthFirstSearchPbcComponents(self, startAtom, parentAtom):
        self.visitedAtoms[startAtom - 1] = 1
        self.atomsInPbcComponent[self.numOfPbcComponents - 1].append(startAtom)
        self.pbcComponentOfAtom[startAtom - 1] = self.numOfPbcComponents - 1
        myNeighbors = self.getNeighbors(startAtom)
        for bondAtom in myNeighbors:
            self.exploreNeighbor(bondAtom, startAtom)

    def getNeighbors(self, atomId):
        return self.bond_atom[atomId]

    def exploreNeighbor(self, bondAtom, startAtom):
        if not self.visitedAtoms[bondAtom - 1]:
            if self.neighborNotAcrossPbcs(bondAtom, startAtom):
                self.depthFirstSearchPbcComponents(bondAtom, startAtom)

    def neighborNotAcrossPbcs(self, bondAtom, startAtom):
        startAtomX = self.getAtomX(startAtom)
        bondAtomX = self.getAtomX(bondAtom)
        for axis in range(3):
            if abs(startAtomX[axis] - bondAtomX[axis]) > self.pbcDistanceThreshold:
                return False
        return True

    def buildGraphAmongPbcComponents(self):
        self.graph = MetaGraph()
        self.edgeList = []
        for pbcComponentId in range(self.numOfPbcComponents):
            self.pbcComponentId = pbcComponentId
            self.connectPbcComponentToItsNeighbors()

    def connectPbcComponentToItsNeighbors(self):
        for atomId in self.atomsInPbcComponent[self.pbcComponentId]:
            atomNeighbors = self.getNeighbors(atomId)
            for atomNeighbor in atomNeighbors:
                self.addEdgeToPbcComponentOfAtomNeighborIfNecessary(
                    atomId, atomNeighbor
                )

    def addEdgeToPbcComponentOfAtomNeighborIfNecessary(self, atomId, atomNeighbor):
        self.neighbors = {}
        neighborPbcComponent = self.pbcComponentOfAtom[atomNeighbor - 1]
        atomIdX = self.getAtomX(atomId)
        neighborX = self.getAtomX(atomNeighbor)
        for axis in range(3):
            boundary = self.findBoundary(atomIdX[axis], neighborX[axis])
            self.updateNeighbors(neighborPbcComponent, axis, boundary)
        self.addEdges()

    def updateNeighbors(self, neighborPbcComponent, axis, boundary):
        if neighborPbcComponent not in self.neighbors:
            self.neighbors[neighborPbcComponent] = []
        self.neighbors[neighborPbcComponent].append(boundary)  # 1, -1 or 0

    def findBoundary(self, atomIdX, neighborX):
        if abs(atomIdX - neighborX) > self.pbcDistanceThreshold:
            if atomIdX < neighborX:
                return -1  # "lo"
            return 1  # "hi"
        return 0  # no pbc crossing in this dimension

    def addEdges(self):
        for neighborPbcComponent in self.neighbors:
            neighborBoundaries = self.neighbors[neighborPbcComponent]
            if neighborBoundaries != [0, 0, 0]:
                self.addEdge(neighborPbcComponent, neighborBoundaries)

    def addEdge(self, neighborPbcComponent, neighborBoundaries):
        edge = (
            self.pbcComponentId,
            neighborPbcComponent,
            neighborBoundaries[0],
            neighborBoundaries[1],
            neighborBoundaries[2],
        )
        if edge not in self.edgeList:
            self.graph.add_edge(edge[0], edge[1], self.createTag(neighborBoundaries))
            self.edgeList.append(edge)
            self.edgeList.append(self.reverseEdge(edge))

    def createTag(self, neighborBoundaries):
        tag = EdgeTag(0, 0, 0)
        for axis, boundary in enumerate(neighborBoundaries):
            if boundary != 0:
                tag[axis] = boundary
        return tag

    def reverseEdge(self, edge):
        return (edge[1], edge[0], -edge[2], -edge[3], -edge[4])

    def findPercolatingMolecules(self):
        self.numPercolatingMolecules = 0
        self.numChains = 0
        self.numSheets = 0
        self.largestMoleculePercolating = False
        self.largest_dimension = 0

        # self.graph.reduce()

        numComponents, dimensionalityList = self.graph.find_stable_loops()

        num_components, component_data = self.graph.get_components()

        # A component, in this context, is not a pbcComponent in the sense
        # meant before. A component is a strongly connected of the metagraph,
        # made of multiple connected nodes. Each node is then a pbcComponent.

        for componentNum in range(numComponents):
            self.largest_dimension = max(
                dimensionalityList[componentNum], self.largest_dimension
            )

            if dimensionalityList[componentNum] == 1:
                self.numChains += 1

            if dimensionalityList[componentNum] == 2:
                self.numSheets += 1

            if dimensionalityList[componentNum] == 3:
                self.numPercolatingMolecules += 1
                component_nodes = component_data[componentNum]
                node_id, node_data = component_nodes[0]  # get node_id of first node

                if not self.largestMoleculePercolating:
                    self.updateLargestMoleculePercolating(node_id)

        """if self.numChains + self.numSheets + self.numPercolatingMolecules > 0:
            self.graph.reduce()
            self.graph.dump(
                "dumps/percolation_frame_{0}.dot".format(self.bondRecord.timeStep)
            )"""

    def updateLargestMoleculePercolating(self, pbcComponentId):
        atomId = self.atomsInPbcComponent[pbcComponentId][0]
        if atomId in self.largestMolecule:
            self.largestMoleculePercolating = True

    # REDUCED MOLECULAR WEIGHT AVERAGE (RMW)

    def computeRMW(self):
        if self.numMolecules == 1:
            self.RMW = 0.0
        else:
            greatestMass = max(self.molecularMasses)
            self.RMW = (self.totalMass - greatestMass) / (self.numMolecules - 1)

    # TWO LARGEST MOLECULAR MASSES OVER TIME

    def getTwoGreatestMolecularMasses(self):
        self.molecularMasses.sort(reverse=True)  # descending order
        self.greatestTwoMasses = [0, 0]
        self.getGreatestMass()
        self.getSecondGreatestMass()

    def getGreatestMass(self):
        self.greatestTwoMasses[0] = self.getMolecularMassFromListIndex(0)

    def getMolecularMassFromListIndex(self, index):
        return self.molecularMasses[index]

    def getSecondGreatestMass(self):
        if self.numMolecules == 1:
            self.greatestTwoMasses[1] = 0.0
        else:
            i = self.getIndexOfSecondGreatestMass()
            self.greatestTwoMasses[1] = self.getMolecularMassFromListIndex(i)

    def getIndexOfSecondGreatestMass(self):
        i = 1
        while self.ithMolecularMassIsWithinNumericalErrorOfMaxMolecularMass(
            self.molecularMasses[i]
        ):
            i += 1
        return i

    def ithMolecularMassIsWithinNumericalErrorOfMaxMolecularMass(self, ithMass):
        max = self.greatestTwoMasses[0]
        epsilon = 0.5
        difference = max - ithMass
        return absDifferenceLessThanEpsilon(difference, epsilon)


def absDifferenceLessThanEpsilon(difference, epsilon):
    if difference < epsilon and difference > (-epsilon):
        return True
    else:
        return False
