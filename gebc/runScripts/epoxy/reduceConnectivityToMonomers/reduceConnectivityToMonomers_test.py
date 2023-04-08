"""Testing network reduction from the full-atom description of the adjancency list to
a coarser, monomer-only description, which retains the correct cross-linking 
connectvity.

The test graph g1 is shown in the pdf file in this folder."""


from os import chdir, path


def reduceBondAdjListToMonomers(bondAdjacencyList, atomTypes):
    """Given the usual bond adjacency list containing all atoms, reduce its 
    granularity by grouping all atom vertices belonging to a single monomer
    into a single vertex"""
    monomersAdjList = {}
    for atomId in bondAdjacencyList:
        if atomTypes[atomId] == "S":
            DFSsearchForNeighboringO(atomId, atomId, atomTypes, monomersAdjList, bondAdjacencyList, {})
    return monomersAdjList


def DFSsearchForNeighboringO(originAtomId, startAtomId, atomTypes, monomersAdjList, bondAdjacencyList, visitedAtoms):
    """Begin DFS from S. Return when encountering either H or O. 
    If O, mark O as neighbor of S and S as neighbor of O."""
    visitedAtoms[startAtomId] = 1
    if atomTypes[startAtomId] == "H":
        return
    if atomTypes[startAtomId] == "O":
        addNeighbor(originAtomId, startAtomId, monomersAdjList)
        addNeighbor(startAtomId, originAtomId, monomersAdjList)
        return
    for childAtomId in bondAdjacencyList[startAtomId]:
        if childAtomId in visitedAtoms:
            continue
        DFSsearchForNeighboringO(originAtomId, childAtomId, atomTypes, monomersAdjList, bondAdjacencyList, visitedAtoms)

    
def addNeighbor(originAtom, newNeighborOfOrigin, monomersAdjList):
    if originAtom in monomersAdjList:
        monomersAdjList[originAtom].append(newNeighborOfOrigin)
    else:
        monomersAdjList[originAtom] = [newNeighborOfOrigin]


def adjListAndAtomTypes(g):
    bondAdjacencyList = {}
    atomTypes = {}
    for atom in g:
        bondAdjacencyList[atom[0]] = g[atom]
        atomTypes[atom[0]] = atom[1]
    return bondAdjacencyList, atomTypes


# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

# test adjacency list reduction
g1 = {}
g1[(0,"C")] = [1]
g1[(1,"O")] = [0, 2]
g1[(2,"N")] = [1, 3, 23]
g1[(3,"O")] = [2, 4]
g1[(4,"N")] = [3, 5, 7]
g1[(5,"O")] = [4, 6]
g1[(6,"C")] = [5]
g1[(7,"S")] = [4, 8]
g1[(8,"N")] = [7, 9, 25]
g1[(9,"O")] = [8, 10]
g1[(10,"N")] = [9, 11, 13]
g1[(11,"O")] = [10, 12]
g1[(12,"C")] = [11]
g1[(13,"S")] = [10, 14]
g1[(14,"N")] = [13, 15, 17]
g1[(15,"O")] = [14, 16]
g1[(16,"C")] = [15]
g1[(17,"O")] = [14, 18]
g1[(18,"N")] = [17, 19, 24]
g1[(19,"S")] = [18, 20]
g1[(20,"N")] = [19, 21, 26]
g1[(21,"O")] = [20, 22]
g1[(22,"N")] = [21, 23, 27]
g1[(23,"S")] = [22, 2]
g1[(24,"H")] = [18]
g1[(25,"H")] = [8]
g1[(26,"H")] = [20]
g1[(27,"H")] = [22]

bondAdjacencyList, atomTypes = adjListAndAtomTypes(g1)
reducedAdjList = reduceBondAdjListToMonomers(bondAdjacencyList, atomTypes)
for vertex in reducedAdjList:
    print(vertex, reducedAdjList[vertex], sep="   ")