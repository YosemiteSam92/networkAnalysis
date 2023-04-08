from geodesicEdgeBetweenneessCentrality\
    import GeodesicEdgeBetweennessCentrality as gebc
from helper_json import loadJson, printJson
from helper_dict import renumberKeysAndValuesFrom0

"""Calculate gebc values for the Voronoi vertices of the cell
   growth simulation. This network was reconstructed by
   reconstructVoronoiConnectivity.py using raw output from
   the simulation. See ../../data/cellGrowth."""

dataPath = "../../data/cellGrowth/simulation3/"
resultsPath = "../../results/cellGrowth/simulation3/"

radius = "24"
adjList_fileName = dataPath + "vertexAdjList_withinRadius"\
    + radius + ".json"
adjList = loadJson(adjList_fileName)

gebcAnalyzer = gebc(adjList, parallel=True, sparseLabels=True)

print("Printing gebc results")
fileName = "voronoiVertices_withinRadius" \
    + radius + "_unnormalized_gebc"

#gebcAnalyzer.printJson(resultsPath + fileName + ".json", normalize=True)
gebcAnalyzer.printJson(resultsPath + fileName + ".json", normalize=False)

print("The end")
