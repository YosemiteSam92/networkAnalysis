import numpy as np
from graph.graph_structs import MetaGraph, EdgeTag, EdgeEntry


def splitRecordEntry(l):
    line = l.split(" ")
    return int(line[1]), int(line[2]), int(line[3])


class MolecularGraphFromMol2:
    def __init__(self, mol2Molecule):
        self.initAttributes(mol2Molecule)
        self.buildMolecularGraph()

    def initAttributes(self, mol2Molecule):
        self.numBonds = mol2Molecule.numBonds
        self.numAtoms = mol2Molecule.numAtoms
        self.bonds = mol2Molecule.triposBond
        self.initBondAtom()
        self.initAtomX(mol2Molecule)  # atomic coordinates

    def initBondAtom(self):
        self.bond_atom = {}  # key: atom id, value: list of atoms bonded to id
        for i in range(1, self.numAtoms + 1):
            self.bond_atom[i] = []

    def initAtomX(self, mol2Molecule):
        self.atomX = {}
        for atom in mol2Molecule.triposAtom:
            id = atom[0]
            self.atomX[id] = self.getAtomCoords(atom)

    def getAtomCoords(self, atom):
        return (atom[1], atom[2], atom[3])

    def buildMolecularGraph(self):
        # vertex = atom; neighbors of vertex = bonding neighbors of atom
        for i in range(self.numBonds):
            id1, id2 = self.getBondPartners(i)
            self.createBond(id1, id2)

    def getBondPartners(self, i):
        return self.bonds[i][0], self.bonds[i][1]

    def createBond(self, id1, id2):
        self.bond_atom[id1].append(id2)
        self.bond_atom[id2].append(id1)

    # PERCOLATING MOLECULE 4
    # find out whether a molecule both spans pbcs in all dimensions.
    # A pbc component containing atom i is composed of all those atoms that can
    # be reached from i via bonds that do not cross any pb.

    def isAnyMoleculePercolating(self):
        self.findAllPbcComponents()
        self.buildGraphAmongPbcComponents()
        self.findPercolatingMolecules()

    def findAllPbcComponents(self):
        self.initAttributesPbcComponents()
        for startAtom in range(1, self.numAtoms + 1):
            self.findPbcComponentContainingStartAtom(startAtom)

    def initAttributesPbcComponents(self):
        self.pbcDistanceThreshold = 0.75
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

    def getAtomX(self, atomId):
        return self.atomX[atomId]

    def buildGraphAmongPbcComponents(self):
        print("self.pbcComponentOfAtom:", self.pbcComponentOfAtom)
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
        print("\n\natomId:", atomId)
        print("neighbor atom:", atomNeighbor)
        self.neighbors = {}
        neighborPbcComponent = self.pbcComponentOfAtom[atomNeighbor - 1]
        atomIdX = self.getAtomX(atomId)
        neighborX = self.getAtomX(atomNeighbor)
        for axis in range(3):
            boundary = self.findBoundary(atomIdX[axis], neighborX[axis])
            self.updateNeighbors(neighborPbcComponent, axis, boundary)

        print("neighbor dict:", self.neighbors)
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
            print("initial node:", self.pbcComponentId)
            print("final node:", neighborPbcComponent)
            print("boundary labels:", neighborBoundaries)
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
        print("\tedge list:", self.edgeList)
        if edge not in self.edgeList:
            print("\tAdding new edge to edge list:", edge)
            self.graph.add_edge(edge[0], edge[1], self.createTag(neighborBoundaries))
            self.edgeList.append(edge)
            self.edgeList.append(self.reverseEdge(edge))

    def createTag(self, neighborBoundaries):
        tag = EdgeTag(0, 0, 0)
        for axis, boundary in enumerate(neighborBoundaries):
            if boundary != 0:
                tag.modify(axis, boundary)
        return tag

    def reverseEdge(self, edge):
        return (edge[1], edge[0], -edge[2], -edge[3], -edge[4])

    def findPercolatingMolecules(self):
        self.numPercolatingMolecules = 0
        self.largestMoleculePercolating = False
        # self.graph.reduce()

        self.numComponents, dimensionalityList = self.graph.find_stable_loops()

        num_components, component_data = self.graph.get_components()

        for componentNum in range(self.numComponents):
            if dimensionalityList[componentNum] == 1:
                self.numPercolatingMolecules += 1

                print("")
                print("Component #{0} data:".format(componentNum))

                component_nodes = component_data[componentNum]
                print("Nodes in component:", component_nodes)

                node_id, node_data = component_nodes[0]
                print("First node in component:", node_id, node_data)

                has_loop = True
                
        print("")
        print("has_loop:", has_loop)
        print("num of components:", self.numComponents)
        print("num of percolating components:", self.numPercolatingMolecules)

        # loop searching but without rotational invariance
        # num_components, has_loop_x = self.graph.find_loops(0)
        # num_components, has_loop_y = self.graph.find_loops(1)
        # num_components, has_loop_z = self.graph.find_loops(2)
        #
        # print("\n\nnum_components = number of molecules:", num_components)
        # print("has_loop_x:", has_loop_x)
        # print("has_loop_y:", has_loop_y)
        # print("has_loop_z:", has_loop_z)

        # for component_num in range(num_components):
        #     if has_loop_x[component_num] and has_loop_y[component_num]: #and has_loop_z[component_num]:
        #         self.numPercolatingMolecules += 1
