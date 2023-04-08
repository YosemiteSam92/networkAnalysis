"""
Given an atom id (possibly involved in bond breakage), determine the monomer
to which it belongs. Monomers are identified by a unique atom representative, as 
stored in the largestMolecularGroups_monomersOnly folder of this repo.

Then, look up the gebc of the monomer.
"""

from lammpsData import LammpsData
from helper_json import loadJson, printJson
from os import chdir, path


def findMonomer(atomId1, atomId2):
    els = (data.getChemicalElement(atomId1),
           data.getChemicalElement(atomId2))
    print(els)
    if els == ("C", "C"):
        # must be DGEBA
        # look for nearest central C
        return (
                DFSsearchForDGEBAcentralC(atomId1, {}, data.bond_atom),
                "DGEBA")


def DFSsearchForDGEBAcentralC(startAtomId, visitedAtoms, bondAdjacencyList):
    """Origin atom is in DGEBA, thus quit the search at 
    either central C or N"""
    visitedAtoms[startAtomId] = 1
    if isDGEBAcentralC(startAtomId):
        return startAtomId
    if data.getChemicalElement(startAtomId) == "N":
        return
    for childAtomId in bondAdjacencyList[startAtomId]:
        if childAtomId in visitedAtoms:
            continue
        centralCid = DFSsearchForDGEBAcentralC(
            childAtomId, visitedAtoms, bondAdjacencyList)
        if centralCid != None:
            return centralCid



def isDGEBAcentralC(atomId):
    if data.getChemicalElement(atomId) == "C":
        if areThere4Cneighbors(data.bond_atom[atomId]):
            return True
    return False


def areThere4Cneighbors(neighborIds):
    """The only DGEBA atom satisfying this condition 
    is the central one."""
    count = 0
    for neighborId in neighborIds:
        if data.getChemicalElement(neighborId) == "C":
            count += 1
    if count == 4:
        return True



# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

# data files
dataFiles = [
    "data.equilibriate_298K_large1_127ns",
    "data.equilibriate_298K_large2_152ns",
    "data.equilibriate_298K_large3_152ns",
    "data.equilibriate_298K_large4_152ns",
    "data.equilibriate_298K_large5_107ns"
]

# adjacency lists with monomers only
monomersAdjLists_fileNames = ["largestGroup_monomers_large"
                              + str(i) for i in range(1,6)]
gebcResults_fileNames = ["largestGroup_monomers_large" + str(i)
                         + "_gebc" for i in range(1,6)]


sysNum = 1
# atom ids of atoms involved in bond breakage
bondAtom1 = 48407
bondAtom2 = 48410

sysNum = 2
bondAtom1 = 87607
bondAtom2 = 87610

sysNum = 3
bondAtom1 = 30473
bondAtom2 = 30476

sysNum = 4
bondAtom1 = 41645
bondAtom2 = 41648

sysNum = 5
bondAtom1 = 16606
bondAtom2 = 16609

dataFile = dataFiles[sysNum - 1]
data = LammpsData("../../lammpsData_annealedEquilibrated/"+dataFile)

# monomersAdjList = loadJson("../../largestMolecularGroups_monomersOnly/"
#                            + monomersAdjLists_fileNames[sysNum - 1] + ".json")

gebc = loadJson("../../results/" + gebcResults_fileNames[sysNum - 1] + ".json")

monomerRepresentativeAtomId, monomerName = findMonomer(bondAtom1, bondAtom2)

monomerRepresentativeAtomId = str(monomerRepresentativeAtomId)
monomerGebc = gebc[monomerRepresentativeAtomId]
print(monomerGebc/max(list(gebc.values())))

