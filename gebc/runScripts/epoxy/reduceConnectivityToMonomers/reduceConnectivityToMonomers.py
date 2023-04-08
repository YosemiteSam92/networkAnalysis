from lammpsData import LammpsData
from os import chdir, path
from sys import setrecursionlimit
from helper_json import printJson
from helper_graph import extractLargestConnectedComponent

"""Given the usual bond adjacency list containing all atoms, reduce its 
    granularity by grouping all atom vertices belonging to a single monomer
    (either DGEBA or DDS) into a single vertex. This operation preserves the 
    cross-linking connectivity of the network, so that no information relevant
    for GEBC calculation is lost (since the topolofy of every monomer 
    is a linear chain). However, the cost of computing GEBC will decrease by
    two orders of magnitude (from about 127000 atoms in the largest molecule
    to about 3000 monomers).
    
    The monomer ids of the reduced network will be given by the atom ids of
    S for a DDS monomer and central C for a DGEBA monomer.
    
    Since the search begins from S atoms and looks for the central C of the
    nearest DGEBA neighbor, unreacted, single DGEBA and DDS monomers will be
    ignored, but as soon as a molecular group contains a cross-link, it will
    be counted. Thus, the largest molecular group still has to be extracted 
    (see below)."""

def reduceBondAdjListToMonomers(bondAdjacencyList, data):
    monomersAdjList = {}
    for atomId in bondAdjacencyList:
        if data.getChemicalElement(atomId) == "S":
            DFSsearchForNeighboringDGEBAcentralC(atomId, atomId, data, monomersAdjList, bondAdjacencyList, {})
    return monomersAdjList


def DFSsearchForNeighboringDGEBAcentralC(originAtomId, startAtomId, data, 
                        monomersAdjList, bondAdjacencyList, visitedAtoms):
    """Begin DFS from S. Return when encountering either H or DGEBA central C. 
    If C, mark C as neighbor of S and S as neighbor of C."""
    visitedAtoms[startAtomId] = 1
    if data.getChemicalElement(startAtomId) == "H":
        return
    if isDGEBAcentralC(startAtomId, data):
        addNeighbor(originAtomId, startAtomId, monomersAdjList)
        addNeighbor(startAtomId, originAtomId, monomersAdjList)
        return
    for childAtomId in bondAdjacencyList[startAtomId]:
        if childAtomId in visitedAtoms:
            continue
        DFSsearchForNeighboringDGEBAcentralC(originAtomId, childAtomId, data, 
                        monomersAdjList, bondAdjacencyList, visitedAtoms)


def isDGEBAcentralC(atomId, data):
    if data.getChemicalElement(atomId) == "C":
        if areThere4Cneighbors(data.bond_atom[atomId], data):
            return True
    return False


def areThere4Cneighbors(neighborIds, data):
    """The only DGEBA atom satisfying this condition 
    is the central one."""
    count = 0
    for neighborId in neighborIds:
        if data.getChemicalElement(neighborId) == "C":
            count += 1
    if count == 4:
        # print("########")
        # for neighborId in neighborIds:
        #     print(data.atomT[neighborId - 1])
        return True


def addNeighbor(originAtom, newNeighborOfOrigin, monomersAdjList):
    if originAtom in monomersAdjList:
        monomersAdjList[originAtom].append(newNeighborOfOrigin)
    else:
        monomersAdjList[originAtom] = [newNeighborOfOrigin]


"""Main program"""
# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

#setrlimit(RLIMIT_STACK, (10 ** 6, -1))
setrecursionlimit(10 ** 6)

dataFiles = [
    "data.equilibriate_298K_large1_127ns",
    "data.equilibriate_298K_large2_152ns",
    "data.equilibriate_298K_large3_152ns",
    "data.equilibriate_298K_large4_152ns",
    "data.equilibriate_298K_large5_107ns"
]

outPath = "../../largestMolecularGroups_monomersOnly/largestGroup_monomers_large"
sysNum = 1

for dataFile in dataFiles:
    data = LammpsData("../../lammpsData_annealedEquilibrated/"+dataFile)
    print("Analizying "+dataFile)
    # data.removeBondsAcrossPbcs()
    largestGroup = extractLargestConnectedComponent(data.bond_atom)
    reducedAdjList = reduceBondAdjListToMonomers(largestGroup, data)
    printJson(reducedAdjList, outPath + str(sysNum) + ".json")
    sysNum += 1


""" 
Results for system 1

Sizes of connected components before reduction:
[126354, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 29, 29]
i.e. the percolating cluster, 12 lone DGEBA molecules and 2 lone DDS molecules. 
Size of largest connected component after reduction:
2986
This is sensible, becuase the previous search had found exactly 
14 unreacted monomers. Since there are initially 3000 monomers in the system,
we would expect 3000-14=2986 monomers in the largest group, as we found.
"""
