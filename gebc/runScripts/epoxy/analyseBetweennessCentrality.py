from geodesicEdgeBetweenneessCentrality\
    import GeodesicEdgeBetweennessCentrality as gebc
from os import chdir, path
from helper_json import loadJson, printJson
from helper_dict import renumberKeysAndValuesFrom0


def printGebc(gebcRaw, conversion, fileName):
    """
    gebcRaw is a 1d numpy array
    apply conversion map to the indices of gebcRaw
    create a dict and print it in json format
    """
    gebc = {}
    for index, value in enumerate(gebcRaw):
        gebc[conversion[index]] = value
    printJson(gebc, fileName)


# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

# adjacency lists of the largest molecular group in the system,
# coarse-grained to include monomers only, which means that only
# information regarding the cross-linking connectivity of the
# network was retained

# The gebc class expects vertices to be labeled incrementally and
# consecutively (this was a short-sighted decision).
# Because the monomer-only adjancency list is no-longer consecutive,
# some renumbering is necessary

adjLists_fileNames = ["largestGroup_monomers_large" + str(i)
                      for i in range(1, 6)]

for fileName in adjLists_fileNames:
    adjList = loadJson("../../data/epoxy/largestMolecularGroups_monomersOnly/"
                       + fileName + ".json")
    print("Analizying "+fileName)
    adjList, consecutiveLabelsToSparseLabels =\
        renumberKeysAndValuesFrom0(adjList)
    gebcAnalyzer = gebc(adjList, parallel=True)

    print("Printing "+fileName)
    # gebcAnalyzer.print("../results/" + fileName + "_gebc_safety.dat")
    # print in json format: monomer id, monomer gebc
    printGebc(gebcAnalyzer.gebc,
              consecutiveLabelsToSparseLabels,
              "../../results/epoxy" + fileName + "_gebc.json")

print("The end")
